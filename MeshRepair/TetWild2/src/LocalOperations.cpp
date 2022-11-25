// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "stdafx.h"
#include "LocalOperations.h"
#include "../external/Predicates.h"
#include "../external/Rational.h"
#include "LocalOperations2.h"
#include <Utils/lowLevel.h>

namespace floatTetWild
{
	const bool use_old_energy = false;
}  // namespace floatTetWild
using floatTetWild::Scalar;

int floatTetWild::get_opp_t_id( const Mesh& mesh, int t_id, int j )
{
	std::vector<int> pair;
	setIntersection( mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 1 ) % 4 ] ].connTets, mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 2 ) % 4 ] ].connTets,
	  mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 3 ) % 4 ] ].connTets, pair );
	if( pair.size() == 2 )
		return pair[ 0 ] == t_id ? pair[ 1 ] : pair[ 0 ];
	return OPP_T_ID_BOUNDARY;
}

void floatTetWild::set_opp_t_id( Mesh& mesh, int t_id, int j )
{
	auto& t = mesh.tets[ t_id ];
	const int jp1 = mod4( j + 1 );
	const int jp2 = mod4( j + 2 );
	const int jp3 = mod4( j + 3 );
	assert( ( j + 1 ) % 4 == jp1 );
	assert( ( j + 2 ) % 4 == jp2 );
	assert( ( j + 3 ) % 4 == jp3 );
	static std::vector<int> pair;
	pair.clear();
	setIntersection( mesh.tet_vertices[ t[ jp1 ] ].connTets, mesh.tet_vertices[ t[ jp2 ] ].connTets, mesh.tet_vertices[ t[ jp3 ] ].connTets, pair );
	if( pair.size() == 2 )
	{
		int opp_t_id = pair[ 0 ] == t_id ? pair[ 1 ] : pair[ 0 ];
		t.opp_t_ids[ j ] = opp_t_id;
		auto& opp_t = mesh.tets[ opp_t_id ];
		for( int k = 0; k < 4; k++ )
		{
			if( opp_t[ k ] != t[ jp1 ] && opp_t[ k ] != t[ jp2 ] && opp_t[ k ] != t[ jp3 ] )
			{
				opp_t.opp_t_ids[ k ] = t_id;
				break;
			}
		}
	}
}

void floatTetWild::get_all_edges( const Mesh& mesh, EdgesSet& edges )
{
	edges.clear();
	edges.reserve( mesh.tets.size() * 6 );

	for( unsigned int i = 0; i < mesh.tets.size(); i++ )
	{
		if( mesh.tets[ i ].is_removed )
			continue;
		for( int j = 0; j < 3; j++ )
		{
			edges.addSorted( mesh.tets[ i ][ 0 ], mesh.tets[ i ][ j + 1 ] );
			edges.addSorted( mesh.tets[ i ][ j + 1 ], mesh.tets[ i ][ mod3( j + 1 ) + 1 ] );
		}
	}
	edges.sortUnique();
}

void floatTetWild::get_all_edges( const Mesh& mesh, const std::vector<int>& t_ids, EdgesSet& edges, bool skip_freezed )
{
	for( unsigned int i = 0; i < t_ids.size(); i++ )
	{
		auto& t = mesh.tets[ t_ids[ i ] ];
		for( int j = 0; j < 3; j++ )
		{
			if( skip_freezed )
			{
				if( !mesh.tet_vertices[ t[ 0 ] ].isFreezed() && !mesh.tet_vertices[ t[ j + 1 ] ].isFreezed() )
					edges.addSorted( t[ 0 ], t[ j + 1 ] );
				if( !mesh.tet_vertices[ t[ j + 1 ] ].isFreezed() && !mesh.tet_vertices[ mod3( j + 1 ) + 1 ].isFreezed() )
					edges.addSorted( t[ j + 1 ], t[ mod3( j + 1 ) + 1 ] );
			}
			else
			{
				edges.addSorted( t[ 0 ], t[ j + 1 ] );
				edges.addSorted( t[ j + 1 ], t[ mod3( j + 1 ) + 1 ] );
			}
		}
	}
	edges.sortUnique();
}

Scalar floatTetWild::get_edge_length( const Mesh& mesh, int v1_id, int v2_id )
{
	const double* const p1 = mesh.tet_vertices[ v1_id ].pos.data();
	const double* const p2 = mesh.tet_vertices[ v2_id ].pos.data();
	// The MeshVertex structure has more fields after the positions, can use unaligned full-vector loads
	const __m256d v1 = _mm256_loadu_pd( p1 );
	const __m256d v2 = _mm256_loadu_pd( p2 );

	const __m256d dist = _mm256_sub_pd( v1, v2 );
	return AvxMath::vector3Length( dist );
}

Scalar floatTetWild::get_edge_length_2( const Mesh& mesh, int v1_id, int v2_id )
{
	const double* const p1 = mesh.tet_vertices[ v1_id ].pos.data();
	const double* const p2 = mesh.tet_vertices[ v2_id ].pos.data();
	// The MeshVertex structure has more fields after the positions, can use unaligned full-vector loads
	const __m256d v1 = _mm256_loadu_pd( p1 );
	const __m256d v2 = _mm256_loadu_pd( p2 );

	const __m256d dist = _mm256_sub_pd( v1, v2 );
	return AvxMath::vector3DotScalar( dist, dist );
}

bool floatTetWild::is_bbox_edge( const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids )
{
	if( !mesh.tet_vertices[ v1_id ].isBoundingBox() || !mesh.tet_vertices[ v2_id ].isBoundingBox() )
		return false;

	for( int t_id : n12_t_ids )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( mesh.tets[ t_id ][ j ] != v1_id && mesh.tets[ t_id ][ j ] != v2_id && mesh.tets[ t_id ].is_bbox_fs[ j ] != NOT_BBOX )
				return true;
		}
	}
	return false;
}

bool floatTetWild::is_surface_edge( const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids )
{
	if( !mesh.tet_vertices[ v1_id ].isSurface() || !mesh.tet_vertices[ v2_id ].isSurface() )
		return false;

	for( int t_id : n12_t_ids )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( mesh.tets[ t_id ][ j ] != v1_id && mesh.tets[ t_id ][ j ] != v2_id && mesh.tets[ t_id ].is_surface_fs[ j ] != NOT_SURFACE )
				return true;
		}
	}
	return false;
}

bool floatTetWild::is_boundary_edge( const Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree )
{
	if( !mesh.tet_vertices[ v1_id ].isBoundary() || !mesh.tet_vertices[ v2_id ].isBoundary() )
		return false;

	std::vector<GEO2::vec3>& ps = mesh.isBoundaryEdgeBuffers.points;
	ps.clear();

	ps.push_back( GEO2::vec3( mesh.tet_vertices[ v1_id ].pos[ 0 ], mesh.tet_vertices[ v1_id ].pos[ 1 ], mesh.tet_vertices[ v1_id ].pos[ 2 ] ) );
	Scalar l = get_edge_length( mesh, v1_id, v2_id );
	const int N = l / mesh.params.dd + 1;
	ps.push_back( GEO2::vec3( mesh.tet_vertices[ v2_id ][ 0 ], mesh.tet_vertices[ v2_id ][ 1 ], mesh.tet_vertices[ v2_id ][ 2 ] ) );

	// TODO: think about reworking this code.
	// Possible to replace AABB tree with BVH, they are often used for ray tracing and can usually query with rays, not just points
	// Here we query these AABB trees with a line segment, pretty much same thing as a ray.
	const __m256d p0 = AvxMath::loadDouble3( ps[ 0 ].data() );
	const __m256d p1 = AvxMath::loadDouble3( ps[ 1 ].data() );

	for( int j = 1; j < N - 1; j++ )
	{
		const double inv = (double)j / (double)N;
		const __m256d pos = AvxMath::lerpFast( p1, p0, inv );
		AvxMath::storeDouble3( ps.emplace_back().data(), pos );
	}

	if( !mesh.is_input_all_inserted )
		return !tree.is_out_tmp_b_envelope( ps, mesh.params.eps_2 );
	else
		return !tree.is_out_b_envelope( ps, mesh.params.eps_2 );
}

bool floatTetWild::is_valid_edge( const Mesh& mesh, int v1_id, int v2_id )
{
	if( mesh.tet_vertices[ v1_id ].isRemoved() || mesh.tet_vertices[ v2_id ].isRemoved() )
		return false;
	return true;
}

bool floatTetWild::is_valid_edge( const Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& n12_t_ids )
{
	if( mesh.tet_vertices[ v1_id ].isRemoved() || mesh.tet_vertices[ v2_id ].isRemoved() )
		return false;
	if( n12_t_ids.empty() )
		return false;
	return true;
}

bool floatTetWild::is_isolate_surface_point( const Mesh& mesh, int v_id )
{
	for( int t_id : mesh.tet_vertices[ v_id ].connTets )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( mesh.tets[ t_id ][ j ] != v_id && mesh.tets[ t_id ].is_surface_fs[ j ] != NOT_SURFACE )
				return false;
		}
	}

	return true;
}

bool floatTetWild::is_point_out_envelope( const Mesh& mesh, const Vector3& p, const AABBWrapper& tree )
{
	GEO2::index_t prev_facet;
	return tree.is_out_sf_envelope( p, mesh.params.eps_2, prev_facet );
}

bool floatTetWild::is_point_out_boundary_envelope( const Mesh& mesh, const Vector3& p, const AABBWrapper& tree )
{
	if( mesh.is_input_all_inserted )
		return false;

	GEO2::index_t prev_facet;
	return tree.is_out_tmp_b_envelope( p, mesh.params.eps_2, prev_facet );
}

namespace
{
	inline void collect2verts( double* rdi, const double* p0, const double* p1 )
	{
		// p0.xy
		__m128d v = _mm_loadu_pd( p0 );
		_mm_storeu_pd( rdi, v );
		// p0.z, p1.x
		v = _mm_load_sd( p0 + 2 );
		v = _mm_loadh_pd( v, p1 );
		_mm_storeu_pd( rdi + 2, v );
		// p1.yz
		v = _mm_loadu_pd( p1 + 1 );
		_mm_storeu_pd( rdi + 4, v );
	}

	inline void collect4verts( std::array<double, 12>& arr, const double* p0, const double* p1, const double* p2, const double* p3 )
	{
		collect2verts( &arr[ 0 ], p0, p1 );
		collect2verts( &arr[ 6 ], p2, p3 );
	}
}  // namespace

Scalar floatTetWild::get_quality( const Mesh& mesh, const MeshTet& t )
{
	std::array<Scalar, 12> T;
	const auto& indices = t.indices;
	collect4verts( T, mesh.tet_vertices[ indices[ 0 ] ].pos.data(), mesh.tet_vertices[ indices[ 1 ] ].pos.data(), mesh.tet_vertices[ indices[ 2 ] ].pos.data(),
	  mesh.tet_vertices[ indices[ 3 ] ].pos.data() );
	return AMIPS_energy( T );
}

Scalar floatTetWild::get_quality( const Mesh& mesh, int t_id )
{
	return get_quality( mesh, mesh.tets[ t_id ] );
}

Scalar floatTetWild::get_quality( const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3 )
{
	std::array<Scalar, 12> T;
	collect4verts( T, v0.data(), v1.data(), v2.data(), v3.data() );
	return AMIPS_energy( T );
}

Scalar floatTetWild::get_quality( const MeshVertex& v0, const MeshVertex& v1, const MeshVertex& v2, const MeshVertex& v3 )
{
	return get_quality( v0.pos, v1.pos, v2.pos, v3.pos );
}

void __declspec( noinline ) floatTetWild::get_max_avg_energy( const Mesh& mesh, Scalar& max_energy, Scalar& avg_energy )
{
	max_energy = 0;
	avg_energy = 0;
	int cnt = 0;
	for( auto& t : mesh.tets )
	{
		if( t.is_removed )
			continue;
		if( t.quality > max_energy )
			max_energy = t.quality;
		avg_energy += t.quality;
		cnt++;
	}
	avg_energy /= cnt;
}

Scalar floatTetWild::get_mid_energy( const Mesh& mesh )
{
	std::vector<Scalar> tmp;
	for( auto& t : mesh.tets )
	{
		if( t.is_removed )
			continue;
		tmp.push_back( t.quality );
	}
	std::sort( tmp.begin(), tmp.end() );
	return tmp[ tmp.size() / 2 ];
}

bool floatTetWild::is_inverted( const Mesh& mesh, int t_id )
{
	const auto& indices = mesh.tets[ t_id ].indices;
	if( Predicates::orient_3d( mesh.tet_vertices[ indices[ 0 ] ].pos, mesh.tet_vertices[ indices[ 1 ] ].pos, mesh.tet_vertices[ indices[ 2 ] ].pos,
		  mesh.tet_vertices[ indices[ 3 ] ].pos ) == eOrientation::Positive )
		return false;
	return true;
}

namespace
{
	// clang-format off
	alignas( 16 ) static const std::array<uint8_t, 16> isInvertedIndices = 
	{
		4, 1, 2, 3, // new_p, v1, v2, v3
		0, 4, 2, 3, // v0, new_p, v2, v3
		0, 1, 4, 3, // v0, v1, new_p, v3
		0, 1, 2, 4	// v0, v1, v2, new_p
	};
	// clang-format on
}  // namespace

bool floatTetWild::is_inverted( const Mesh& mesh, int t_id, int j, const Vector3& new_p )
{
	assert( j >= 0 && j < 4 );

	const Vector4i& tet = mesh.tets[ t_id ].indices;
	std::array<const Vector3*, 5> sourcePointers = {
	  &mesh.tet_vertices[ tet[ 0 ] ].pos, &mesh.tet_vertices[ tet[ 1 ] ].pos, &mesh.tet_vertices[ tet[ 2 ] ].pos, &mesh.tet_vertices[ tet[ 3 ] ].pos, &new_p };

	const uint8_t* lookup = isInvertedIndices.data() + 4 * j;
	const Vector3** rsi = sourcePointers.data();
	eOrientation ori = Predicates::orient_3d( *rsi[ lookup[ 0 ] ], *rsi[ lookup[ 1 ] ], *rsi[ lookup[ 2 ] ], *rsi[ lookup[ 3 ] ] );
	return ori != eOrientation::Positive;
}

bool floatTetWild::is_inverted( const MeshVertex& v0, const MeshVertex& v1, const MeshVertex& v2, const MeshVertex& v3 )
{
	if( Predicates::orient_3d( v0.pos, v1.pos, v2.pos, v3.pos ) == eOrientation::Positive )
		return false;
	return true;
}

bool floatTetWild::is_inverted( const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3 )
{
	if( Predicates::orient_3d( v0, v1, v2, v3 ) == eOrientation::Positive )
		return false;
	return true;
}

bool floatTetWild::is_degenerate( const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& v3 )
{
	if( Predicates::orient_3d( v0, v1, v2, v3 ) == eOrientation::Zero )
		return true;
	return false;
}

bool floatTetWild::is_out_boundary_envelope( const Mesh& mesh, int v_id, const Vector3& new_pos, const AABBWrapper& tree )
{
	if( mesh.is_input_all_inserted )
		return false;
	if( !mesh.tet_vertices[ v_id ].isCut() )
		return false;

	GEO2::index_t prev_facet;
	if( tree.is_out_tmp_b_envelope( new_pos, mesh.params.eps_2 / 100, prev_facet ) )
		return true;

	std::vector<int> tmp_b_v_ids;
	for( int t_id : mesh.tet_vertices[ v_id ].connTets )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( mesh.tets[ t_id ][ j ] != v_id && mesh.tets[ t_id ].is_surface_fs[ j ] <= 0 )
			{
				for( int k = 0; k < 3; k++ )
				{
					int b_v_id = mesh.tets[ t_id ][ ( j + 1 + k ) % 4 ];
					if( b_v_id != v_id && mesh.tet_vertices[ b_v_id ].isBoundary() )
						tmp_b_v_ids.push_back( b_v_id );
				}
			}
		}
	}
	vector_unique( tmp_b_v_ids );

	std::vector<int> b_v_ids;
	b_v_ids.reserve( tmp_b_v_ids.size() );
	for( int b_v_id : tmp_b_v_ids )
	{
		if( is_boundary_edge( mesh, v_id, b_v_id, tree ) )
			b_v_ids.push_back( b_v_id );
	}
	if( b_v_ids.empty() )
		return false;

	std::vector<GEO2::vec3> ps;
	ps.push_back( GEO2::vec3( new_pos[ 0 ], new_pos[ 1 ], new_pos[ 2 ] ) );
	int p0_id = 0;
	for( int b_v_id : b_v_ids )
	{
		Scalar l = get_edge_length( mesh, v_id, b_v_id );
		int N = l / mesh.params.dd + 1;
		ps.push_back( GEO2::vec3( mesh.tet_vertices[ b_v_id ][ 0 ], mesh.tet_vertices[ b_v_id ][ 1 ], mesh.tet_vertices[ b_v_id ][ 2 ] ) );
		int p1_id = ps.size() - 1;
		for( Scalar j = 1; j < N - 1; j++ )
		{
			ps.push_back( ps[ p0_id ] * ( j / N ) + ps[ p1_id ] * ( 1 - j / N ) );
		}
	}
	return tree.is_out_tmp_b_envelope( ps, mesh.params.eps_2 / 100, prev_facet );
}

bool floatTetWild::is_out_envelope( Mesh& mesh, int v_id, const Vector3& new_pos, const AABBWrapper& tree )
{
	GEO2::index_t prev_facet;
	if( tree.is_out_sf_envelope( new_pos, mesh.params.eps_2, prev_facet ) )
		return true;

	for( int t_id : mesh.tet_vertices[ v_id ].connTets )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( mesh.tets[ t_id ][ j ] != v_id && mesh.tets[ t_id ].is_surface_fs[ j ] <= 0 )
			{
				std::array<Vector3, 3> vs;
				for( int k = 0; k < 3; k++ )
				{
					if( mesh.tets[ t_id ][ mod4( j + 1 + k ) ] == v_id )
						vs[ k ] = new_pos;
					else
						vs[ k ] = mesh.tet_vertices[ mesh.tets[ t_id ][ mod4( j + 1 + k ) ] ].pos;
				}
				bool is_out = sample_triangle_and_check_is_out( vs, mesh.params.dd, mesh.params.eps_2, tree, prev_facet );
				if( is_out )
					return true;
			}
		}
	}
	return false;
}

void floatTetWild::sample_triangle( const std::array<Vector3, 3>& vs, std::vector<GEO2::vec3>& ps, Scalar sampling_dist )
{
	Scalar sqrt3_2 = std::sqrt( 3 ) / 2;

	std::array<Scalar, 3> ls;
	for( int i = 0; i < 3; i++ )
	{
		ls[ i ] = ( vs[ i ] - vs[ mod3( i + 1 ) ] ).squaredNorm();
	}
	auto min_max = std::minmax_element( ls.begin(), ls.end() );
	int min_i = min_max.first - ls.begin();
	int max_i = min_max.second - ls.begin();
	Scalar N = sqrt( ls[ max_i ] ) / sampling_dist;
	if( N <= 1 )
	{
		for( int i = 0; i < 3; i++ )
			ps.push_back( GEO2::vec3( vs[ i ][ 0 ], vs[ i ][ 1 ], vs[ i ][ 2 ] ) );
		return;
	}
	if( N == int( N ) )
		N -= 1;

	GEO2::vec3 v0( vs[ max_i ][ 0 ], vs[ max_i ][ 1 ], vs[ max_i ][ 2 ] );
	GEO2::vec3 v1( vs[ mod3( max_i + 1 ) ][ 0 ], vs[ mod3( max_i + 1 ) ][ 1 ], vs[ mod3( max_i + 1 ) ][ 2 ] );
	GEO2::vec3 v2( vs[ mod3( max_i + 2 ) ][ 0 ], vs[ mod3( max_i + 2 ) ][ 1 ], vs[ mod3( max_i + 2 ) ][ 2 ] );

	GEO2::vec3 n_v0v1 = GEO2::normalize( v1 - v0 );
	for( int n = 0; n <= N; n++ )
	{
		ps.push_back( v0 + n_v0v1 * sampling_dist * n );
	}
	ps.push_back( v1 );

	Scalar h = GEO2::distance( GEO2::dot( ( v2 - v0 ), ( v1 - v0 ) ) * ( v1 - v0 ) / ls[ max_i ] + v0, v2 );
	int M = h / ( sqrt3_2 * sampling_dist );
	if( M < 1 )
	{
		ps.push_back( v2 );
		return;
	}

	GEO2::vec3 n_v0v2 = GEO2::normalize( v2 - v0 );
	GEO2::vec3 n_v1v2 = GEO2::normalize( v2 - v1 );
	Scalar tan_v0, tan_v1, sin_v0, sin_v1;
	sin_v0 = GEO2::length( GEO2::cross( ( v2 - v0 ), ( v1 - v0 ) ) ) / ( GEO2::distance( v0, v2 ) * GEO2::distance( v0, v1 ) );
	tan_v0 = GEO2::length( GEO2::cross( ( v2 - v0 ), ( v1 - v0 ) ) ) / GEO2::dot( ( v2 - v0 ), ( v1 - v0 ) );
	tan_v1 = GEO2::length( GEO2::cross( ( v2 - v1 ), ( v0 - v1 ) ) ) / GEO2::dot( ( v2 - v1 ), ( v0 - v1 ) );
	sin_v1 = GEO2::length( GEO2::cross( ( v2 - v1 ), ( v0 - v1 ) ) ) / ( GEO2::distance( v1, v2 ) * GEO2::distance( v0, v1 ) );

	for( int m = 1; m <= M; m++ )
	{
		int n = sqrt3_2 / tan_v0 * m + 0.5;
		int n1 = sqrt3_2 / tan_v0 * m;
		if( m % 2 == 0 && n == n1 )
		{
			n += 1;
		}
		GEO2::vec3 v0_m = v0 + m * sqrt3_2 * sampling_dist / sin_v0 * n_v0v2;
		GEO2::vec3 v1_m = v1 + m * sqrt3_2 * sampling_dist / sin_v1 * n_v1v2;
		if( GEO2::distance( v0_m, v1_m ) <= sampling_dist )
			break;

		Scalar delta_d = ( ( n + ( m % 2 ) / 2.0 ) - m * sqrt3_2 / tan_v0 ) * sampling_dist;
		GEO2::vec3 v = v0_m + delta_d * n_v0v1;
		int N1 = GEO2::distance( v, v1_m ) / sampling_dist;
		//        ps.push_back(v0_m);
		for( int i = 0; i <= N1; i++ )
		{
			ps.push_back( v + i * n_v0v1 * sampling_dist );
		}
		//        ps.push_back(v1_m);
	}
	ps.push_back( v2 );

	// sample edges
	N = sqrt( ls[ mod3( max_i + 1 ) ] ) / sampling_dist;
	if( N > 1 )
	{
		if( N == int( N ) )
			N -= 1;
		GEO2::vec3 n_v1v2 = GEO2::normalize( v2 - v1 );
		for( int n = 1; n <= N; n++ )
		{
			ps.push_back( v1 + n_v1v2 * sampling_dist * n );
		}
	}

	N = sqrt( ls[ mod3( max_i + 2 ) ] ) / sampling_dist;
	if( N > 1 )
	{
		if( N == int( N ) )
			N -= 1;
		GEO2::vec3 n_v2v0 = GEO2::normalize( v0 - v2 );
		for( int n = 1; n <= N; n++ )
		{
			ps.push_back( v2 + n_v2v0 * sampling_dist * n );
		}
	}
}

bool floatTetWild::sample_triangle_and_check_is_out(
  const std::array<Vector3, 3>& vs, Scalar sampling_dist, Scalar eps_2, const AABBWrapper& tree, GEO2::index_t& prev_facet )
{
	GEO2::vec3 nearest_point;
	double sq_dist = std::numeric_limits<double>::max();

	Scalar sqrt3_2 = std::sqrt( 3 ) / 2;

	std::array<Scalar, 3> ls;
	for( int i = 0; i < 3; i++ )
	{
		ls[ i ] = ( vs[ i ] - vs[ mod3( i + 1 ) ] ).squaredNorm();
	}
	auto min_max = std::minmax_element( ls.begin(), ls.end() );
	int min_i = min_max.first - ls.begin();
	int max_i = min_max.second - ls.begin();
	Scalar N = sqrt( ls[ max_i ] ) / sampling_dist;
	if( N <= 1 )
	{
		for( int i = 0; i < 3; i++ )
		{
			if( tree.is_out_sf_envelope( vs[ i ], eps_2, prev_facet, sq_dist, nearest_point ) )
				return true;
		}
		return false;
	}
	if( N == int( N ) )
		N -= 1;

	GEO2::vec3 v0( vs[ max_i ][ 0 ], vs[ max_i ][ 1 ], vs[ max_i ][ 2 ] );
	GEO2::vec3 v1( vs[ mod3( max_i + 1 ) ][ 0 ], vs[ mod3( max_i + 1 ) ][ 1 ], vs[ mod3( max_i + 1 ) ][ 2 ] );
	GEO2::vec3 v2( vs[ mod3( max_i + 2 ) ][ 0 ], vs[ mod3( max_i + 2 ) ][ 1 ], vs[ mod3( max_i + 2 ) ][ 2 ] );

	GEO2::vec3 n_v0v1 = GEO2::normalize( v1 - v0 );
	for( int n = 0; n <= N; n++ )
	{
		if( tree.is_out_sf_envelope( v0 + n_v0v1 * sampling_dist * n, eps_2, prev_facet, sq_dist, nearest_point ) )
			return true;
	}
	if( tree.is_out_sf_envelope( v1, eps_2, prev_facet, sq_dist, nearest_point ) )
		return true;

	Scalar h = GEO2::distance( GEO2::dot( ( v2 - v0 ), ( v1 - v0 ) ) * ( v1 - v0 ) / ls[ max_i ] + v0, v2 );
	int M = h / ( sqrt3_2 * sampling_dist );
	if( M < 1 )
		return tree.is_out_sf_envelope( v2, eps_2, prev_facet, sq_dist, nearest_point );

	GEO2::vec3 n_v0v2 = GEO2::normalize( v2 - v0 );
	GEO2::vec3 n_v1v2 = GEO2::normalize( v2 - v1 );
	Scalar tan_v0, tan_v1, sin_v0, sin_v1;
	sin_v0 = GEO2::length( GEO2::cross( ( v2 - v0 ), ( v1 - v0 ) ) ) / ( GEO2::distance( v0, v2 ) * GEO2::distance( v0, v1 ) );
	tan_v0 = GEO2::length( GEO2::cross( ( v2 - v0 ), ( v1 - v0 ) ) ) / GEO2::dot( ( v2 - v0 ), ( v1 - v0 ) );
	tan_v1 = GEO2::length( GEO2::cross( ( v2 - v1 ), ( v0 - v1 ) ) ) / GEO2::dot( ( v2 - v1 ), ( v0 - v1 ) );
	sin_v1 = GEO2::length( GEO2::cross( ( v2 - v1 ), ( v0 - v1 ) ) ) / ( GEO2::distance( v1, v2 ) * GEO2::distance( v0, v1 ) );

	for( int m = 1; m <= M; m++ )
	{
		int n = sqrt3_2 / tan_v0 * m + 0.5;
		int n1 = sqrt3_2 / tan_v0 * m;
		if( m % 2 == 0 && n == n1 )
			n++;

		GEO2::vec3 v0_m = v0 + m * sqrt3_2 * sampling_dist / sin_v0 * n_v0v2;
		GEO2::vec3 v1_m = v1 + m * sqrt3_2 * sampling_dist / sin_v1 * n_v1v2;
		if( GEO2::distance( v0_m, v1_m ) <= sampling_dist )
			break;

		Scalar delta_d = ( ( n + ( m % 2 ) / 2.0 ) - m * sqrt3_2 / tan_v0 ) * sampling_dist;
		GEO2::vec3 v = v0_m + delta_d * n_v0v1;
		int N1 = GEO2::distance( v, v1_m ) / sampling_dist;
		for( int i = 0; i <= N1; i++ )
		{
			if( tree.is_out_sf_envelope( v + i * n_v0v1 * sampling_dist, eps_2, prev_facet, sq_dist, nearest_point ) )
				return true;
		}
	}

	if( tree.is_out_sf_envelope( v2, eps_2, prev_facet, sq_dist, nearest_point ) )
		return true;

	// sample edges
	N = sqrt( ls[ mod3( max_i + 1 ) ] ) / sampling_dist;
	if( N > 1 )
	{
		if( N == int( N ) )
			N -= 1;
		GEO2::vec3 n_v1v2 = GEO2::normalize( v2 - v1 );
		for( int n = 1; n <= N; n++ )
		{
			if( tree.is_out_sf_envelope( v1 + n_v1v2 * sampling_dist * n, eps_2, prev_facet, sq_dist, nearest_point ) )
				return true;
		}
	}

	N = sqrt( ls[ mod3( max_i + 2 ) ] ) / sampling_dist;
	if( N > 1 )
	{
		if( N == int( N ) )
			N -= 1;
		GEO2::vec3 n_v2v0 = GEO2::normalize( v0 - v2 );
		for( int n = 1; n <= N; n++ )
		{
			if( tree.is_out_sf_envelope( v2 + n_v2v0 * sampling_dist * n, eps_2, prev_facet, sq_dist, nearest_point ) )
				return true;
		}
	}

	return false;
}

void floatTetWild::get_new_tet_slots( Mesh& mesh, int n, std::vector<int>& new_conn_tets )
{
	int cnt = 0;
	for( int i = mesh.t_empty_start; i < mesh.tets.size(); i++ )
	{
		if( mesh.tets[ i ].is_removed )
		{
			new_conn_tets.push_back( i );
			cnt++;
			if( cnt == n )
			{
				mesh.t_empty_start = i + 1;
				break;
			}
		}
	}
	if( cnt < n )
	{
		for( int i = 0; i < n - cnt; i++ )
			new_conn_tets.push_back( mesh.tets.size() + i );
		mesh.tets.resize( mesh.tets.size() + n - cnt );
		mesh.t_empty_start = mesh.tets.size();
	}
}

void floatTetWild::set_intersection( const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::vector<int>& v )
{
	if( s2.size() < s1.size() )
	{
		set_intersection( s2, s1, v );
		return;
	}
	v.clear();
	v.reserve( std::min( s1.size(), s2.size() ) );
	for( int x : s1 )
	{
		if( s2.count( x ) )
		{
			v.push_back( x );
		}
	}
}

void floatTetWild::set_intersection( const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, std::unordered_set<int>& v )
{
	if( s2.size() < s1.size() )
	{
		set_intersection( s2, s1, v );
		return;
	}
	v.clear();
	v.reserve( std::min( s1.size(), s2.size() ) );
	for( int x : s1 )
	{
		if( s2.count( x ) )
		{
			v.insert( x );
		}
	}
}

void floatTetWild::set_intersection(
  const std::unordered_set<int>& s1, const std::unordered_set<int>& s2, const std::unordered_set<int>& s3, std::vector<int>& v )
{
	if( s2.size() < s1.size() && s2.size() < s1.size() )
	{
		set_intersection( s2, s1, s3, v );
		return;
	}

	if( s3.size() < s1.size() && s3.size() < s2.size() )
	{
		set_intersection( s3, s1, s2, v );
		return;
	}

	assert( s1.size() <= s2.size() );
	assert( s1.size() <= s3.size() );

	v.clear();
	v.reserve( s1.size() );
	for( int x : s1 )
	{
		if( s2.count( x ) && s3.count( x ) )
		{
			v.push_back( x );
			if( v.size() == 2 )
				break;
		}
	}
}

void floatTetWild::set_intersection( const std::vector<int>& s11, const std::vector<int>& s22, std::vector<int>& v )
{
	std::vector<int> s1 = s11;
	std::vector<int> s2 = s22;
	std::sort( s1.begin(), s1.end() );
	std::sort( s2.begin(), s2.end() );
	std::set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter( v ) );
}

void floatTetWild::set_intersection( const std::vector<int>& s11, const std::vector<int>& s22, const std::vector<int>& s33, std::vector<int>& v )
{
	std::vector<int> s1 = s11;
	std::vector<int> s2 = s22;
	std::vector<int> s3 = s33;
	std::sort( s1.begin(), s1.end() );
	std::sort( s2.begin(), s2.end() );
	std::sort( s3.begin(), s3.end() );
	std::set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter( v ) );
	auto it = std::set_intersection( v.begin(), v.end(), s3.begin(), s3.end(), v.begin() );
	v.resize( it - v.begin() );
}

void floatTetWild::set_intersection_sorted( const std::vector<int>& s1, const std::vector<int>& s2, const std::vector<int>& s3, std::vector<int>& v )
{
	std::set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter( v ) );
	auto it = std::set_intersection( v.begin(), v.end(), s3.begin(), s3.end(), v.begin() );
	v.resize( it - v.begin() );
}

void floatTetWild::pausee( std::string msg )
{
	return;
}

namespace
{
	inline double pow2( double x )
	{
		return x * x;
	}

	inline double cubicRoot( double x )
	{
		return floatTetWild::cbrt( x );
	}
}  // namespace

Scalar floatTetWild::AMIPS_energy( const std::array<Scalar, 12>& T )
{
	Scalar res = AMIPS_energy_aux( T );
	if( use_old_energy )
	{
		return res;
	}

	if( res > 1e8 )
	{
		if( is_degenerate(
			  Vector3( T[ 0 ], T[ 1 ], T[ 2 ] ), Vector3( T[ 3 ], T[ 4 ], T[ 5 ] ), Vector3( T[ 6 ], T[ 7 ], T[ 8 ] ), Vector3( T[ 9 ], T[ 10 ], T[ 11 ] ) ) )
		{
			pausee( "energy computation degenerate found!!!" );
			return std::numeric_limits<double>::infinity();
		}

		std::array<triwild::Rational, 12> r_T;
		for( int j = 0; j < 12; j++ )
			r_T[ j ] = T[ j ];
		const triwild::Rational twothird = triwild::Rational( 2 ) / triwild::Rational( 3 );
		triwild::Rational tmp =
		  ( ( -r_T[ 1 + 2 ] + r_T[ 1 + 5 ] ) * r_T[ 1 + 1 ] + r_T[ 1 + 2 ] * r_T[ 1 + 7 ] + ( r_T[ 1 + -1 ] - r_T[ 1 + 5 ] ) * r_T[ 1 + 4 ] -
			r_T[ 1 + -1 ] * r_T[ 1 + 7 ] ) *
			r_T[ 1 + 9 ] +
		  ( ( r_T[ 1 + 2 ] - r_T[ 1 + 5 ] ) * r_T[ 1 + 0 ] - r_T[ 1 + 2 ] * r_T[ 1 + 6 ] + ( -r_T[ 1 + -1 ] + r_T[ 1 + 5 ] ) * r_T[ 1 + 3 ] +
			r_T[ 1 + -1 ] * r_T[ 1 + 6 ] ) *
			r_T[ 1 + 10 ] +
		  ( -r_T[ 1 + 2 ] * r_T[ 1 + 7 ] + ( -r_T[ 1 + 8 ] + r_T[ 1 + 5 ] ) * r_T[ 1 + 4 ] + r_T[ 1 + 8 ] * r_T[ 1 + 7 ] ) * r_T[ 1 + 0 ] +
		  ( r_T[ 1 + 2 ] * r_T[ 1 + 6 ] + ( r_T[ 1 + 8 ] - r_T[ 1 + 5 ] ) * r_T[ 1 + 3 ] - r_T[ 1 + 8 ] * r_T[ 1 + 6 ] ) * r_T[ 1 + 1 ] +
		  ( r_T[ 1 + 3 ] * r_T[ 1 + 7 ] - r_T[ 1 + 4 ] * r_T[ 1 + 6 ] ) * ( r_T[ 1 + -1 ] - r_T[ 1 + 8 ] );
		if( tmp == 0 )
			return std::numeric_limits<double>::infinity();

		auto res_r =
		  triwild::Rational( 27 ) / 16 * pow( tmp, -2 ) *
		  pow( r_T[ 1 + 9 ] * r_T[ 1 + 9 ] + ( -twothird * r_T[ 1 + 0 ] - twothird * r_T[ 1 + 3 ] - twothird * r_T[ 1 + 6 ] ) * r_T[ 1 + 9 ] +
				 r_T[ 1 + 10 ] * r_T[ 1 + 10 ] + ( -twothird * r_T[ 1 + 1 ] - twothird * r_T[ 1 + 4 ] - twothird * r_T[ 1 + 7 ] ) * r_T[ 1 + 10 ] +
				 r_T[ 1 + 0 ] * r_T[ 1 + 0 ] + ( -twothird * r_T[ 1 + 3 ] - twothird * r_T[ 1 + 6 ] ) * r_T[ 1 + 0 ] + r_T[ 1 + 1 ] * r_T[ 1 + 1 ] +
				 ( -twothird * r_T[ 1 + 4 ] - twothird * r_T[ 1 + 7 ] ) * r_T[ 1 + 1 ] + r_T[ 1 + 2 ] * r_T[ 1 + 2 ] +
				 ( -twothird * r_T[ 1 + -1 ] - twothird * r_T[ 1 + 8 ] - twothird * r_T[ 1 + 5 ] ) * r_T[ 1 + 2 ] + r_T[ 1 + 3 ] * r_T[ 1 + 3 ] -
				 twothird * r_T[ 1 + 3 ] * r_T[ 1 + 6 ] + r_T[ 1 + 4 ] * r_T[ 1 + 4 ] - twothird * r_T[ 1 + 4 ] * r_T[ 1 + 7 ] + r_T[ 1 + 5 ] * r_T[ 1 + 5 ] +
				 ( -twothird * r_T[ 1 + -1 ] - twothird * r_T[ 1 + 8 ] ) * r_T[ 1 + 5 ] - twothird * r_T[ 1 + -1 ] * r_T[ 1 + 8 ] +
				 r_T[ 1 + -1 ] * r_T[ 1 + -1 ] + r_T[ 1 + 8 ] * r_T[ 1 + 8 ] + r_T[ 1 + 6 ] * r_T[ 1 + 6 ] + r_T[ 1 + 7 ] * r_T[ 1 + 7 ],
			3 );
		return cubicRoot( res_r.to_double() );
	}
	else
		return res;
}

static Scalar AMIPS_energy_aux_v1( const std::array<Scalar, 12>& T )
{
	Scalar helper_0[ 12 ];
	helper_0[ 0 ] = T[ 0 ];
	helper_0[ 1 ] = T[ 1 ];
	helper_0[ 2 ] = T[ 2 ];
	helper_0[ 3 ] = T[ 3 ];
	helper_0[ 4 ] = T[ 4 ];
	helper_0[ 5 ] = T[ 5 ];
	helper_0[ 6 ] = T[ 6 ];
	helper_0[ 7 ] = T[ 7 ];
	helper_0[ 8 ] = T[ 8 ];
	helper_0[ 9 ] = T[ 9 ];
	helper_0[ 10 ] = T[ 10 ];
	helper_0[ 11 ] = T[ 11 ];
	Scalar helper_1 = helper_0[ 2 ];
	Scalar helper_2 = helper_0[ 11 ];
	Scalar helper_3 = helper_0[ 0 ];
	Scalar helper_4 = helper_0[ 3 ];
	Scalar helper_5 = helper_0[ 9 ];
	Scalar helper_6 = 0.577350269189626 * helper_3 - 1.15470053837925 * helper_4 + 0.577350269189626 * helper_5;
	Scalar helper_7 = helper_0[ 1 ];
	Scalar helper_8 = helper_0[ 4 ];
	Scalar helper_9 = helper_0[ 7 ];
	Scalar helper_10 = helper_0[ 10 ];
	Scalar helper_11 = 0.408248290463863 * helper_10 + 0.408248290463863 * helper_7 + 0.408248290463863 * helper_8 - 1.22474487139159 * helper_9;
	Scalar helper_12 = 0.577350269189626 * helper_10 + 0.577350269189626 * helper_7 - 1.15470053837925 * helper_8;
	Scalar helper_13 = helper_0[ 6 ];
	Scalar helper_14 = -1.22474487139159 * helper_13 + 0.408248290463863 * helper_3 + 0.408248290463863 * helper_4 + 0.408248290463863 * helper_5;
	Scalar helper_15 = helper_0[ 5 ];
	Scalar helper_16 = helper_0[ 8 ];
	Scalar helper_17 = 0.408248290463863 * helper_1 + 0.408248290463863 * helper_15 - 1.22474487139159 * helper_16 + 0.408248290463863 * helper_2;
	Scalar helper_18 = 0.577350269189626 * helper_1 - 1.15470053837925 * helper_15 + 0.577350269189626 * helper_2;
	Scalar helper_19 = 0.5 * helper_13 + 0.5 * helper_4;
	Scalar helper_20 = 0.5 * helper_8 + 0.5 * helper_9;
	Scalar helper_21 = 0.5 * helper_15 + 0.5 * helper_16;
	Scalar helper_22 = ( helper_1 - helper_2 ) * ( helper_11 * helper_6 - helper_12 * helper_14 ) -
					   ( -helper_10 + helper_7 ) * ( -helper_14 * helper_18 + helper_17 * helper_6 ) +
					   ( helper_3 - helper_5 ) * ( -helper_11 * helper_18 + helper_12 * helper_17 );
	Scalar res =
	  -( helper_1 * ( -1.5 * helper_1 + 0.5 * helper_2 + helper_21 ) + helper_10 * ( -1.5 * helper_10 + helper_20 + 0.5 * helper_7 ) +
		 helper_13 * ( -1.5 * helper_13 + 0.5 * helper_3 + 0.5 * helper_4 + 0.5 * helper_5 ) +
		 helper_15 * ( 0.5 * helper_1 - 1.5 * helper_15 + 0.5 * helper_16 + 0.5 * helper_2 ) +
		 helper_16 * ( 0.5 * helper_1 + 0.5 * helper_15 - 1.5 * helper_16 + 0.5 * helper_2 ) + helper_2 * ( 0.5 * helper_1 - 1.5 * helper_2 + helper_21 ) +
		 helper_3 * ( helper_19 - 1.5 * helper_3 + 0.5 * helper_5 ) + helper_4 * ( 0.5 * helper_13 + 0.5 * helper_3 - 1.5 * helper_4 + 0.5 * helper_5 ) +
		 helper_5 * ( helper_19 + 0.5 * helper_3 - 1.5 * helper_5 ) + helper_7 * ( 0.5 * helper_10 + helper_20 - 1.5 * helper_7 ) +
		 helper_8 * ( 0.5 * helper_10 + 0.5 * helper_7 - 1.5 * helper_8 + 0.5 * helper_9 ) +
		 helper_9 * ( 0.5 * helper_10 + 0.5 * helper_7 + 0.5 * helper_8 - 1.5 * helper_9 ) ) /
	  cubicRoot( helper_22 * helper_22 );
	return res;
}

namespace
{
	bool areEqualRel( double a, double b, double epsilon )
	{
		const bool f1 = std::isfinite( a );
		const bool f2 = std::isfinite( b );
		if( f1 && f2 )
			return ( std::abs( a - b ) <= epsilon * std::max( std::abs( a ), std::abs( b ) ) );
		return f1 == f2;
	}
}  // namespace

Scalar floatTetWild::AMIPS_energy_aux( const std::array<Scalar, 12>& T )
{
#if 1
	return AMIPS_energy_aux_v4( T );
#else
	const double v1 = AMIPS_energy_aux_v1( T );
	const double v3 = AMIPS_energy_aux_v3( T );
	const double v2 = AMIPS_energy_aux_v4( T );
	// if( v1 != v2 ) __debugbreak(); return v2;
	constexpr double eps = 1E-4;
	if( !areEqualRel( v1, v2, eps ) )
	{
		if( v2 <= 1e8 )
			__debugbreak();
		// When it exceeds 1E+8, the calling code in AMIPS_energy() will use the rational numbers anyway
	}
	return v2;
#endif
}

void floatTetWild::AMIPS_jacobian( const std::array<Scalar, 12>& T, Vector3& result_0 )
{
	Scalar helper_0[ 12 ];
	helper_0[ 0 ] = T[ 0 ];
	helper_0[ 1 ] = T[ 1 ];
	helper_0[ 2 ] = T[ 2 ];
	helper_0[ 3 ] = T[ 3 ];
	helper_0[ 4 ] = T[ 4 ];
	helper_0[ 5 ] = T[ 5 ];
	helper_0[ 6 ] = T[ 6 ];
	helper_0[ 7 ] = T[ 7 ];
	helper_0[ 8 ] = T[ 8 ];
	helper_0[ 9 ] = T[ 9 ];
	helper_0[ 10 ] = T[ 10 ];
	helper_0[ 11 ] = T[ 11 ];
	Scalar helper_1 = helper_0[ 1 ];
	Scalar helper_2 = helper_0[ 10 ];
	Scalar helper_3 = helper_1 - helper_2;
	Scalar helper_4 = helper_0[ 0 ];
	Scalar helper_5 = helper_0[ 3 ];
	Scalar helper_6 = helper_0[ 9 ];
	Scalar helper_7 = 0.577350269189626 * helper_4 - 1.15470053837925 * helper_5 + 0.577350269189626 * helper_6;
	Scalar helper_8 = helper_0[ 2 ];
	Scalar helper_9 = 0.408248290463863 * helper_8;
	Scalar helper_10 = helper_0[ 5 ];
	Scalar helper_11 = 0.408248290463863 * helper_10;
	Scalar helper_12 = helper_0[ 8 ];
	Scalar helper_13 = 1.22474487139159 * helper_12;
	Scalar helper_14 = helper_0[ 11 ];
	Scalar helper_15 = 0.408248290463863 * helper_14;
	Scalar helper_16 = helper_11 - helper_13 + helper_15 + helper_9;
	Scalar helper_17 = 0.577350269189626 * helper_8;
	Scalar helper_18 = 1.15470053837925 * helper_10;
	Scalar helper_19 = 0.577350269189626 * helper_14;
	Scalar helper_20 = helper_17 - helper_18 + helper_19;
	Scalar helper_21 = helper_0[ 6 ];
	Scalar helper_22 = -1.22474487139159 * helper_21 + 0.408248290463863 * helper_4 + 0.408248290463863 * helper_5 + 0.408248290463863 * helper_6;
	Scalar helper_23 = helper_16 * helper_7 - helper_20 * helper_22;
	Scalar helper_24 = -helper_14 + helper_8;
	Scalar helper_25 = 0.408248290463863 * helper_1;
	Scalar helper_26 = helper_0[ 4 ];
	Scalar helper_27 = 0.408248290463863 * helper_26;
	Scalar helper_28 = helper_0[ 7 ];
	Scalar helper_29 = 1.22474487139159 * helper_28;
	Scalar helper_30 = 0.408248290463863 * helper_2;
	Scalar helper_31 = helper_25 + helper_27 - helper_29 + helper_30;
	Scalar helper_32 = helper_31 * helper_7;
	Scalar helper_33 = 0.577350269189626 * helper_1;
	Scalar helper_34 = 1.15470053837925 * helper_26;
	Scalar helper_35 = 0.577350269189626 * helper_2;
	Scalar helper_36 = helper_33 - helper_34 + helper_35;
	Scalar helper_37 = helper_22 * helper_36;
	Scalar helper_38 = helper_4 - helper_6;
	Scalar helper_39 = helper_23 * helper_3 - helper_24 * ( helper_32 - helper_37 ) - helper_38 * ( helper_16 * helper_36 - helper_20 * helper_31 );
	Scalar helper_40 = 1.0 / cubicRoot( pow2( helper_39 ) );
	Scalar helper_41 = 0.707106781186548 * helper_10 - 0.707106781186548 * helper_12;
	Scalar helper_42 = 0.707106781186548 * helper_26 - 0.707106781186548 * helper_28;
	Scalar helper_43 = 0.5 * helper_21 + 0.5 * helper_5;
	Scalar helper_44 = 0.5 * helper_26 + 0.5 * helper_28;
	Scalar helper_45 = 0.5 * helper_10 + 0.5 * helper_12;
	Scalar helper_46 =
	  0.666666666666667 *
	  ( helper_1 * ( -1.5 * helper_1 + 0.5 * helper_2 + helper_44 ) + helper_10 * ( -1.5 * helper_10 + 0.5 * helper_12 + 0.5 * helper_14 + 0.5 * helper_8 ) +
		helper_12 * ( 0.5 * helper_10 - 1.5 * helper_12 + 0.5 * helper_14 + 0.5 * helper_8 ) + helper_14 * ( -1.5 * helper_14 + helper_45 + 0.5 * helper_8 ) +
		helper_2 * ( 0.5 * helper_1 - 1.5 * helper_2 + helper_44 ) + helper_21 * ( -1.5 * helper_21 + 0.5 * helper_4 + 0.5 * helper_5 + 0.5 * helper_6 ) +
		helper_26 * ( 0.5 * helper_1 + 0.5 * helper_2 - 1.5 * helper_26 + 0.5 * helper_28 ) +
		helper_28 * ( 0.5 * helper_1 + 0.5 * helper_2 + 0.5 * helper_26 - 1.5 * helper_28 ) + helper_4 * ( -1.5 * helper_4 + helper_43 + 0.5 * helper_6 ) +
		helper_5 * ( 0.5 * helper_21 + 0.5 * helper_4 - 1.5 * helper_5 + 0.5 * helper_6 ) + helper_6 * ( 0.5 * helper_4 + helper_43 - 1.5 * helper_6 ) +
		helper_8 * ( 0.5 * helper_14 + helper_45 - 1.5 * helper_8 ) ) /
	  helper_39;
	Scalar helper_47 = -0.707106781186548 * helper_21 + 0.707106781186548 * helper_5;
	result_0[ 0 ] = -helper_40 * ( 1.0 * helper_21 - 3.0 * helper_4 +
								   helper_46 * ( helper_41 * ( -helper_1 + helper_2 ) - helper_42 * ( helper_14 - helper_8 ) -
												 ( -helper_17 + helper_18 - helper_19 ) * ( -helper_25 - helper_27 + helper_29 - helper_30 ) +
												 ( -helper_33 + helper_34 - helper_35 ) * ( -helper_11 + helper_13 - helper_15 - helper_9 ) ) +
								   1.0 * helper_5 + 1.0 * helper_6 );
	result_0[ 1 ] = helper_40 * ( 3.0 * helper_1 - 1.0 * helper_2 - 1.0 * helper_26 - 1.0 * helper_28 +
								  helper_46 * ( helper_23 + helper_24 * helper_47 - helper_38 * helper_41 ) );
	result_0[ 2 ] = helper_40 * ( -1.0 * helper_10 - 1.0 * helper_12 - 1.0 * helper_14 +
								  helper_46 * ( -helper_3 * helper_47 - helper_32 + helper_37 + helper_38 * helper_42 ) + 3.0 * helper_8 );
}

void floatTetWild::AMIPS_hessian( const std::array<Scalar, 12>& T, Matrix3& result_0 )
{
	Scalar helper_0[ 12 ];
	helper_0[ 0 ] = T[ 0 ];
	helper_0[ 1 ] = T[ 1 ];
	helper_0[ 2 ] = T[ 2 ];
	helper_0[ 3 ] = T[ 3 ];
	helper_0[ 4 ] = T[ 4 ];
	helper_0[ 5 ] = T[ 5 ];
	helper_0[ 6 ] = T[ 6 ];
	helper_0[ 7 ] = T[ 7 ];
	helper_0[ 8 ] = T[ 8 ];
	helper_0[ 9 ] = T[ 9 ];
	helper_0[ 10 ] = T[ 10 ];
	helper_0[ 11 ] = T[ 11 ];
	Scalar helper_1 = helper_0[ 2 ];
	Scalar helper_2 = helper_0[ 11 ];
	Scalar helper_3 = helper_1 - helper_2;
	Scalar helper_4 = helper_0[ 0 ];
	Scalar helper_5 = 0.577350269189626 * helper_4;
	Scalar helper_6 = helper_0[ 3 ];
	Scalar helper_7 = 1.15470053837925 * helper_6;
	Scalar helper_8 = helper_0[ 9 ];
	Scalar helper_9 = 0.577350269189626 * helper_8;
	Scalar helper_10 = helper_5 - helper_7 + helper_9;
	Scalar helper_11 = helper_0[ 1 ];
	Scalar helper_12 = 0.408248290463863 * helper_11;
	Scalar helper_13 = helper_0[ 4 ];
	Scalar helper_14 = 0.408248290463863 * helper_13;
	Scalar helper_15 = helper_0[ 7 ];
	Scalar helper_16 = 1.22474487139159 * helper_15;
	Scalar helper_17 = helper_0[ 10 ];
	Scalar helper_18 = 0.408248290463863 * helper_17;
	Scalar helper_19 = helper_12 + helper_14 - helper_16 + helper_18;
	Scalar helper_20 = helper_10 * helper_19;
	Scalar helper_21 = 0.577350269189626 * helper_11;
	Scalar helper_22 = 1.15470053837925 * helper_13;
	Scalar helper_23 = 0.577350269189626 * helper_17;
	Scalar helper_24 = helper_21 - helper_22 + helper_23;
	Scalar helper_25 = 0.408248290463863 * helper_4;
	Scalar helper_26 = 0.408248290463863 * helper_6;
	Scalar helper_27 = helper_0[ 6 ];
	Scalar helper_28 = 1.22474487139159 * helper_27;
	Scalar helper_29 = 0.408248290463863 * helper_8;
	Scalar helper_30 = helper_25 + helper_26 - helper_28 + helper_29;
	Scalar helper_31 = helper_24 * helper_30;
	Scalar helper_32 = helper_3 * ( helper_20 - helper_31 );
	Scalar helper_33 = helper_4 - helper_8;
	Scalar helper_34 = 0.408248290463863 * helper_1;
	Scalar helper_35 = helper_0[ 5 ];
	Scalar helper_36 = 0.408248290463863 * helper_35;
	Scalar helper_37 = helper_0[ 8 ];
	Scalar helper_38 = 1.22474487139159 * helper_37;
	Scalar helper_39 = 0.408248290463863 * helper_2;
	Scalar helper_40 = helper_34 + helper_36 - helper_38 + helper_39;
	Scalar helper_41 = helper_24 * helper_40;
	Scalar helper_42 = 0.577350269189626 * helper_1;
	Scalar helper_43 = 1.15470053837925 * helper_35;
	Scalar helper_44 = 0.577350269189626 * helper_2;
	Scalar helper_45 = helper_42 - helper_43 + helper_44;
	Scalar helper_46 = helper_19 * helper_45;
	Scalar helper_47 = helper_41 - helper_46;
	Scalar helper_48 = helper_33 * helper_47;
	Scalar helper_49 = helper_11 - helper_17;
	Scalar helper_50 = helper_10 * helper_40;
	Scalar helper_51 = helper_30 * helper_45;
	Scalar helper_52 = helper_50 - helper_51;
	Scalar helper_53 = helper_49 * helper_52;
	Scalar helper_54 = helper_32 + helper_48 - helper_53;
	Scalar helper_55 = pow2( helper_54 );
	Scalar helper_56 = 1.0 / cubicRoot( helper_55 );
	Scalar helper_57 = 1.0 * helper_27 - 3.0 * helper_4 + 1.0 * helper_6 + 1.0 * helper_8;
	Scalar helper_58 = 0.707106781186548 * helper_13;
	Scalar helper_59 = 0.707106781186548 * helper_15;
	Scalar helper_60 = helper_58 - helper_59;
	Scalar helper_61 = helper_3 * helper_60;
	Scalar helper_62 = 0.707106781186548 * helper_35 - 0.707106781186548 * helper_37;
	Scalar helper_63 = helper_49 * helper_62;
	Scalar helper_64 = helper_47 + helper_61 - helper_63;
	Scalar helper_65 = 1.33333333333333 / helper_54;
	Scalar helper_66 = 1.0 / helper_55;
	Scalar helper_67 = 0.5 * helper_27 + 0.5 * helper_6;
	Scalar helper_68 = -1.5 * helper_4 + helper_67 + 0.5 * helper_8;
	Scalar helper_69 = 0.5 * helper_4 + helper_67 - 1.5 * helper_8;
	Scalar helper_70 = -1.5 * helper_27 + 0.5 * helper_4 + 0.5 * helper_6 + 0.5 * helper_8;
	Scalar helper_71 = 0.5 * helper_27 + 0.5 * helper_4 - 1.5 * helper_6 + 0.5 * helper_8;
	Scalar helper_72 = 0.5 * helper_13 + 0.5 * helper_15;
	Scalar helper_73 = -1.5 * helper_11 + 0.5 * helper_17 + helper_72;
	Scalar helper_74 = 0.5 * helper_11 - 1.5 * helper_17 + helper_72;
	Scalar helper_75 = 0.5 * helper_11 + 0.5 * helper_13 - 1.5 * helper_15 + 0.5 * helper_17;
	Scalar helper_76 = 0.5 * helper_11 - 1.5 * helper_13 + 0.5 * helper_15 + 0.5 * helper_17;
	Scalar helper_77 = 0.5 * helper_35 + 0.5 * helper_37;
	Scalar helper_78 = -1.5 * helper_1 + 0.5 * helper_2 + helper_77;
	Scalar helper_79 = 0.5 * helper_1 - 1.5 * helper_2 + helper_77;
	Scalar helper_80 = 0.5 * helper_1 + 0.5 * helper_2 + 0.5 * helper_35 - 1.5 * helper_37;
	Scalar helper_81 = 0.5 * helper_1 + 0.5 * helper_2 - 1.5 * helper_35 + 0.5 * helper_37;
	Scalar helper_82 = helper_1 * helper_78 + helper_11 * helper_73 + helper_13 * helper_76 + helper_15 * helper_75 + helper_17 * helper_74 +
					   helper_2 * helper_79 + helper_27 * helper_70 + helper_35 * helper_81 + helper_37 * helper_80 + helper_4 * helper_68 +
					   helper_6 * helper_71 + helper_69 * helper_8;
	Scalar helper_83 = 0.444444444444444 * helper_66 * helper_82;
	Scalar helper_84 = helper_66 * helper_82;
	Scalar helper_85 = -helper_32 - helper_48 + helper_53;
	Scalar helper_86 = 1.0 / helper_85;
	Scalar helper_87 = helper_86 / cubicRoot( pow2( helper_85 ) );
	Scalar helper_88 = 0.707106781186548 * helper_6;
	Scalar helper_89 = 0.707106781186548 * helper_27;
	Scalar helper_90 = helper_88 - helper_89;
	Scalar helper_91 = 0.666666666666667 * helper_10 * helper_40 + 0.666666666666667 * helper_3 * helper_90 - 0.666666666666667 * helper_30 * helper_45 -
					   0.666666666666667 * helper_33 * helper_62;
	Scalar helper_92 = -3.0 * helper_11 + 1.0 * helper_13 + 1.0 * helper_15 + 1.0 * helper_17;
	Scalar helper_93 = -helper_11 + helper_17;
	Scalar helper_94 = -helper_1 + helper_2;
	Scalar helper_95 = -helper_21 + helper_22 - helper_23;
	Scalar helper_96 = -helper_34 - helper_36 + helper_38 - helper_39;
	Scalar helper_97 = -helper_42 + helper_43 - helper_44;
	Scalar helper_98 = -helper_12 - helper_14 + helper_16 - helper_18;
	Scalar helper_99 = -0.666666666666667 * helper_60 * helper_94 + 0.666666666666667 * helper_62 * helper_93 + 0.666666666666667 * helper_95 * helper_96 -
					   0.666666666666667 * helper_97 * helper_98;
	Scalar helper_100 = helper_3 * helper_90;
	Scalar helper_101 = helper_33 * helper_62;
	Scalar helper_102 = helper_100 - helper_101 + helper_52;
	Scalar helper_103 = -helper_60 * helper_94 + helper_62 * helper_93 + helper_95 * helper_96 - helper_97 * helper_98;
	Scalar helper_104 = 0.444444444444444 * helper_102 * helper_103 * helper_82 * helper_86 + helper_57 * helper_91 - helper_92 * helper_99;
	Scalar helper_105 =
	  1.85037170770859e-17 * helper_1 * helper_78 + 1.85037170770859e-17 * helper_11 * helper_73 + 1.85037170770859e-17 * helper_13 * helper_76 +
	  1.85037170770859e-17 * helper_15 * helper_75 + 1.85037170770859e-17 * helper_17 * helper_74 + 1.85037170770859e-17 * helper_2 * helper_79 +
	  1.85037170770859e-17 * helper_27 * helper_70 + 1.85037170770859e-17 * helper_35 * helper_81 + 1.85037170770859e-17 * helper_37 * helper_80 +
	  1.85037170770859e-17 * helper_4 * helper_68 + 1.85037170770859e-17 * helper_6 * helper_71 + 1.85037170770859e-17 * helper_69 * helper_8;
	Scalar helper_106 = helper_64 * helper_82 * helper_86;
	Scalar helper_107 = -0.666666666666667 * helper_10 * helper_19 + 0.666666666666667 * helper_24 * helper_30 + 0.666666666666667 * helper_33 * helper_60 -
						0.666666666666667 * helper_49 * helper_90;
	Scalar helper_108 = -3.0 * helper_1 + 1.0 * helper_2 + 1.0 * helper_35 + 1.0 * helper_37;
	Scalar helper_109 = -helper_20 + helper_31 + helper_33 * helper_60 - helper_49 * helper_90;
	Scalar helper_110 = 0.444444444444444 * helper_109 * helper_82 * helper_86;
	Scalar helper_111 = helper_103 * helper_110 + helper_107 * helper_57 - helper_108 * helper_99;
	Scalar helper_112 = -helper_4 + helper_8;
	Scalar helper_113 = -helper_88 + helper_89;
	Scalar helper_114 = -helper_5 + helper_7 - helper_9;
	Scalar helper_115 = -helper_25 - helper_26 + helper_28 - helper_29;
	Scalar helper_116 = helper_82 * helper_86 * ( helper_112 * helper_62 + helper_113 * helper_94 + helper_114 * helper_96 - helper_115 * helper_97 );
	Scalar helper_117 = -helper_100 + helper_101 - helper_50 + helper_51;
	Scalar helper_118 = -helper_102 * helper_110 + helper_107 * helper_92 + helper_108 * helper_91;
	Scalar helper_119 =
	  helper_82 * helper_86 * ( helper_112 * ( -helper_58 + helper_59 ) - helper_113 * helper_93 - helper_114 * helper_98 + helper_115 * helper_95 );
	result_0( 0, 0 ) = helper_56 * ( helper_57 * helper_64 * helper_65 - pow2( helper_64 ) * helper_83 +
									 0.666666666666667 * helper_64 * helper_84 * ( -helper_41 + helper_46 - helper_61 + helper_63 ) + 3.0 );
	result_0( 0, 1 ) = helper_87 * ( helper_104 - helper_105 * helper_35 + helper_106 * helper_91 );
	result_0( 0, 2 ) = helper_87 * ( helper_106 * helper_107 + helper_111 );
	result_0( 1, 0 ) = helper_87 * ( helper_104 + helper_116 * helper_99 );
	result_0( 1, 1 ) = helper_56 * ( -pow2( helper_117 ) * helper_83 + helper_117 * helper_65 * helper_92 + helper_117 * helper_84 * helper_91 + 3.0 );
	result_0( 1, 2 ) = helper_87 * ( -helper_105 * helper_6 - helper_107 * helper_116 + helper_118 );
	result_0( 2, 0 ) = helper_87 * ( -helper_105 * helper_13 + helper_111 + helper_119 * helper_99 );
	result_0( 2, 1 ) = helper_87 * ( helper_118 - helper_119 * helper_91 );
	result_0( 2, 2 ) = helper_56 * ( -helper_108 * helper_109 * helper_65 - 1.11111111111111 * pow2( helper_109 ) * helper_84 + 3.0 );
}

void floatTetWild::vector_unique( std::vector<int>& vec )
{
	std::sort( vec.begin(), vec.end() );
	vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );
}

namespace
{
	// A faster comparison predicate to sort std::array<int, 2>
	// Only delivers equal result to the built-in `operator <` when the integers aren't negative. Good enough for the use case.
	struct Int2Cmp
	{
		inline bool operator()( uint64_t a, uint64_t b )
		{
			// Swap low and high 32-bit pieces, with the rotate instructions
			a = _rotr64( a, 32 );
			b = _rotr64( b, 32 );
			// Now the comparison is identical to the std::lexicographical_compare of the arrays
			return a < b;
		}
	};
}  // namespace

void floatTetWild::vector_unique( std::vector<std::array<int, 2>>& vec )
{
	if( vec.empty() )
		return;

	static_assert( sizeof( std::array<int, 2> ) == sizeof( uint64_t ) );
	// std::vector<std::array<int, 2>> cpy = vec;
	// std::sort( cpy.begin(), cpy.end() );

	// Sort the vector using that faster comparison predicate
	uint64_t* const begin = (uint64_t*)vec.data();
	uint64_t* const end = begin + vec.size();
	std::sort( begin, end, Int2Cmp {} );

	// if( cpy != vec ) __debugbreak();

	// Run std::unique algorithm, pretending these pairs are uint64_t values
	// The default comparison and assignment already do the right thing for these values
	uint64_t* const unique = std::unique( begin, end );

	// Finally, shrink the vector
	vec.erase( vec.begin() + ( unique - begin ), vec.end() );
}

void floatTetWild::vector_unique( std::vector<std::array<int, 3>>& vec )
{
	std::sort( vec.begin(), vec.end() );
	vec.erase( std::unique( vec.begin(), vec.end() ), vec.end() );
}