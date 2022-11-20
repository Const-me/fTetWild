// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "stdafx.h"
#include "FloatTetDelaunay.h"
#include <igl/Timer.h>
#include <iterator>
#include <algorithm>
#include <bitset>
#include "../external/Predicates.h"
#include "LocalOperations.h"
#include "MeshImprovement.h"
#include "../../GeogramDelaunay/GeogramDelaunay.h"
#include "BoolVector.h"

namespace floatTetWild
{
	namespace
	{
		void get_bb_corners( const Parameters& params, const std::vector<Vector3>& vertices, Vector3& min, Vector3& max )
		{
			min = vertices.front();
			max = vertices.front();

			for( size_t j = 0; j < vertices.size(); j++ )
			{
				for( int i = 0; i < 3; i++ )
				{
					min( i ) = std::min( min( i ), vertices[ j ]( i ) );
					max( i ) = std::max( max( i ), vertices[ j ]( i ) );
				}
			}

			//            const Scalar dis = std::max((max - min).minCoeff() * params.box_scale, params.eps_input * 2);
			const Scalar dis = std::max( params.ideal_edge_length, params.eps_input * 2 );
			for( int j = 0; j < 3; j++ )
			{
				min[ j ] -= dis;
				max[ j ] += dis;
			}

			params.logger.logDebug( "min = %g %g %g", min[ 0 ], min[ 1 ], min[ 2 ] );
			params.logger.logDebug( "max = %g %g %g", max[ 0 ], max[ 1 ], max[ 2 ] );
		}

		bool comp( const std::array<int, 4>& a, const std::array<int, 4>& b )
		{
			return std::tuple<int, int, int>( a[ 0 ], a[ 1 ], a[ 2 ] ) < std::tuple<int, int, int>( b[ 0 ], b[ 1 ], b[ 2 ] );
		}

		void match_surface_fs( Mesh& mesh, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, BoolVector& is_face_inserted )
		{
			std::vector<std::array<int, 4>> input_fs( input_faces.size() );
			for( int i = 0; i < input_faces.size(); i++ )
			{
				input_fs[ i ] = { { input_faces[ i ][ 0 ], input_faces[ i ][ 1 ], input_faces[ i ][ 2 ], i } };
				std::sort( input_fs[ i ].begin(), input_fs[ i ].begin() + 3 );
			}
			std::sort( input_fs.begin(), input_fs.end(), comp );

			for( auto& t : mesh.tets )
			{
				for( int j = 0; j < 4; j++ )
				{
					std::array<int, 3> f = { { t[ ( j + 1 ) % 4 ], t[ ( j + 2 ) % 4 ], t[ ( j + 3 ) % 4 ] } };
					std::sort( f.begin(), f.end() );
					auto bounds = std::equal_range( input_fs.begin(), input_fs.end(), std::array<int, 4>( { { f[ 0 ], f[ 1 ], f[ 2 ], -1 } } ), comp );
					for( auto it = bounds.first; it != bounds.second; ++it )
					{
						int f_id = ( *it )[ 3 ];
						is_face_inserted.set( f_id );
					}
				}
			}
		}

		void match_bbox_fs( Mesh& mesh, const Vector3& min, const Vector3& max )
		{
			auto get_bbox_fs = [ & ]( const MeshTet& t, int j )
			{
				std::array<int, 6> cnts = { { 0, 0, 0, 0, 0, 0 } };
				for( int k = 0; k < 3; k++ )
				{
					Vector3& pos = mesh.tet_vertices[ t[ ( j + k + 1 ) % 4 ] ].pos;
					for( int n = 0; n < 3; n++ )
					{
						if( pos[ n ] == min[ n ] )
							cnts[ n * 2 ]++;
						else if( pos[ n ] == max[ n ] )
							cnts[ n * 2 + 1 ]++;
					}
				}
				for( int i = 0; i < cnts.size(); i++ )
				{
					if( cnts[ i ] == 3 )
						return (int8_t)i;
				}
				return NOT_BBOX;
			};

			for( auto& t : mesh.tets )
			{
				for( int j = 0; j < 4; j++ )
				{
					t.is_bbox_fs[ j ] = get_bbox_fs( t, j );
				}
			}
		}

		namespace
		{
			inline double lerpFast( double x, double y, double s )
			{
				return x * ( 1.0 - s ) + y * s;
			}
		}  // namespace

		__m128i computeVoxelPoints(
		  const Vector3& minScalar, const Vector3& maxScalar, const Parameters& params, const AABBWrapper& tree, std::vector<Vector3>& voxels )
		{
			using namespace AvxMath;
			const __m256d min = loadDouble3( minScalar.data() );
			const __m256d max = loadDouble3( maxScalar.data() );
			const __m256d diag = _mm256_sub_pd( max, min );

			__m128i nVoxels = _mm256_cvtpd_epi32( _mm256_div_pd( diag, _mm256_set1_pd( params.bbox_diag_length * params.box_scale ) ) );
			nVoxels = _mm_max_epi32( nVoxels, _mm_set1_epi32( 1 ) );

			const size_t vz = (uint32_t)_mm_extract_epi32( nVoxels, 2 );
			const size_t vy = (uint32_t)_mm_extract_epi32( nVoxels, 1 );
			const size_t vx = (uint32_t)_mm_cvtsi128_si32( nVoxels );

			voxels.clear();
			voxels.resize( ( vx + 1 ) * ( vy + 1 ) * ( vz + 1 ) );
			Vector3* rdi = voxels.data();

			const double sq_distg = 100 * params.eps_2;

			Vector3 voxelDiv;
			storeDouble3( voxelDiv.data(), _mm256_cvtepi32_pd( nVoxels ) );
			Vector3 pos;

			for( size_t z = 0; z <= vz; z++ )
			{
				pos[ 2 ] = lerpFast( minScalar[ 2 ], maxScalar[ 2 ], (double)(int)z / voxelDiv[ 2 ] );

				for( size_t y = 0; y <= vy; y++ )
				{
					pos[ 1 ] = lerpFast( minScalar[ 1 ], maxScalar[ 1 ], (double)(int)y / voxelDiv[ 1 ] );
					for( size_t x = 0; x <= vx; x++ )
					{
						pos[ 0 ] = lerpFast( minScalar[ 0 ], maxScalar[ 0 ], (double)(int)x / voxelDiv[ 0 ] );
						if( tree.get_sq_dist_to_sf( pos ) > sq_distg )
						{
							*rdi = pos;
							rdi++;
						}
					}
				}
			}
			voxels.resize( rdi - voxels.data() );
			return nVoxels;
		}

		__forceinline __m256d loadVertex( const double* vb, size_t vbSize, int si )
		{
			const size_t ui = (size_t)(uint32_t)si;
			const double* const p = vb + ui * 3;
			if( ui + 1 < vbSize )
				return _mm256_loadu_pd( p );  //< This branch is extremely likely to be taken, and saves quite a few instructions
			else
				return AvxMath::loadDouble3( p );
		}

		__forceinline void updateMaxTetraSize( __m256d& acc, const double* vb, size_t vbSize, __m128i elt )
		{
			const __m256d v0 = loadVertex( vb, vbSize, _mm_cvtsi128_si32( elt ) );
			const __m256d v1 = loadVertex( vb, vbSize, _mm_extract_epi32( elt, 1 ) );
			const __m256d v2 = loadVertex( vb, vbSize, _mm_extract_epi32( elt, 2 ) );
			const __m256d v3 = loadVertex( vb, vbSize, _mm_extract_epi32( elt, 3 ) );

			__m256d i = _mm256_min_pd( _mm256_min_pd( v0, v1 ), _mm256_min_pd( v2, v3 ) );
			__m256d ax = _mm256_max_pd( _mm256_max_pd( v0, v1 ), _mm256_max_pd( v2, v3 ) );
			__m256d sz = _mm256_sub_pd( ax, i );
			acc = _mm256_max_pd( acc, sz );
		}
	}  // namespace

	void FloatTetDelaunay::tetrahedralize(
	  const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const AABBWrapper& tree, Mesh& mesh, BoolVector& is_face_inserted )
	{
		const Parameters& params = mesh.params;

		is_face_inserted.resize( input_faces.size(), false );

		Vector3 min, max;
		get_bb_corners( params, input_vertices, min, max );
		mesh.params.bbox_min = min;
		mesh.params.bbox_max = max;

		std::vector<Vector3> voxel_points;
#if PARALLEL_TRIANGLES_INSERTION
		const __m128i voxelCount = computeVoxelPoints( min, max, params, tree, voxel_points );
#else
		computeVoxelPoints( min, max, params, tree, voxel_points );
#endif
		const size_t n_pts = input_vertices.size() + voxel_points.size();

		std::vector<Vector3> allPoints;
		allPoints.resize( n_pts );
		std::copy( input_vertices.begin(), input_vertices.end(), allPoints.begin() );
		std::copy( voxel_points.begin(), voxel_points.end(), allPoints.begin() + input_vertices.size() );

		auto& tet_vertices = mesh.tet_vertices;
		tet_vertices.resize( n_pts );
		for( size_t i = 0; i < n_pts; i++ )
			tet_vertices[ i ].pos = allPoints[ i ];

		// The Delaunay is implemented in Geogram; we consume the algorithm through the usable wrapper exposed by GeogramDelaunay static library
		constexpr bool multithreadedDelaunay = false;
		auto delaunay = iDelaunay::create( multithreadedDelaunay );
		delaunay->compute( n_pts, (const double*)allPoints.data() );

		const size_t countElements = delaunay->countElements();
		auto& tets = mesh.tets;
		tets.resize( countElements );
		const __m128i* const tet2v = delaunay->getElements();
		for( size_t i = 0; i < countElements; i++ )
		{
			__m128i tetra = tet2v[ i ];
			tetra = _mm_shuffle_epi32( tetra, _MM_SHUFFLE( 1, 2, 3, 0 ) );
			_mm_storeu_si128( (__m128i*)&tets[ i ].indices, tetra );

			tet_vertices[ _mm_cvtsi128_si32( tetra ) ].connTets.add( (int)i );
			tet_vertices[ _mm_extract_epi32( tetra, 1 ) ].connTets.add( (int)i );
			tet_vertices[ _mm_extract_epi32( tetra, 2 ) ].connTets.add( (int)i );
			tet_vertices[ _mm_extract_epi32( tetra, 3 ) ].connTets.add( (int)i );
		}

#if PARALLEL_TRIANGLES_INSERTION
		const double* const vertexPointer = (const double*)allPoints.data();
		const size_t vertexCount = tet_vertices.size();
		__m256d maxElementSizeVec = _mm256_setzero_pd();

		for( size_t i = 0; i < countElements; i++ )
		{
			// The above loop calls vector.push_back()
			// The compiler won't inline that, and due to the ABI convention this evicts some of the registers to memory
			// updateMaxTetraSize uses quite a few of these vector registers.
			// It's probably more efficient to run another loop over the same elements, this time without calling any functions from the loop
			updateMaxTetraSize( maxElementSizeVec, vertexPointer, vertexCount, tet2v[ i ] );
		}

		// We gonna divide vectors by that vector; set W lane to 1.0 to avoid division by 0 in the unused lane
		// On some CPUs, these singularity shenanigans have a small performance cost
		maxElementSizeVec = _mm256_blend_pd( maxElementSizeVec, _mm256_set1_pd( 1 ), 0b1000 );
		// Store in the mesh
		mesh.maxTetraSize = maxElementSizeVec;
#endif

		for( int i = 0; i < mesh.tets.size(); i++ )
		{
			if( is_inverted( mesh, i ) )
				throw std::logic_error( "EXIT_INV" );
		}

		// match faces: should be integer with sign
		// match bbox 8 facets: should be -1 and 0~5
		//        match_surface_fs(mesh, input_vertices, input_faces, is_face_inserted);
		match_bbox_fs( mesh, min, max );
	}
}  // namespace floatTetWild
