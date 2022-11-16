/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */

#include "stdafx.h"
#include "mesh_AABB.h"
#include "../Utils/BoundingBox.hpp"
#ifdef __AVX__
#include <Utils/AvxMath.h>
#endif
#include <omp.h>

namespace
{
	using namespace GEO2;
	using namespace floatTetWild;

	/**
	 * \brief Computes the axis-aligned bounding box of a mesh facet.
	 * \param[in] M the mesh
	 * \param[out] B the bounding box of the facet
	 * \param[in] f the index of the facet in mesh \p M
	 */
	void get_facet_bbox( const Mesh& M, Box& B, index_t f )
	{
		const vec3 *p1, *p2, *p3;
		M.getTriangleVertices( f, &p1, &p2, &p3 );

		BoundingBox bb( &p1->x );
		bb.extend( &p2->x );
		bb.extend( &p3->x );
		bb.store( B );
	}

	/**
	 * \brief Computes the maximum node index in a subtree
	 * \param[in] node_index node index of the root of the subtree
	 * \param[in] b first facet index in the subtree
	 * \param[in] e one position past the last facet index in the subtree
	 * \return the maximum node index in the subtree rooted at \p node_index
	 */
	index_t max_node_index( index_t node_index, index_t b, index_t e )
	{
		assert( e > b );
		if( b + 1 == e )
		{
			return node_index;
		}
		index_t m = b + ( e - b ) / 2;
		index_t childl = 2 * node_index;
		index_t childr = 2 * node_index + 1;
		return std::max( max_node_index( childl, b, m ), max_node_index( childr, m, e ) );
	}

	/**
	 * \brief Computes the hierarchy of bounding boxes recursively.
	 * \details This function is generic and can be used to compute
	 *  a bbox hierarchy of arbitrary elements.
	 * \param[in] M the mesh
	 * \param[in] bboxes the array of bounding boxes
	 * \param[in] node_index the index of the root of the subtree
	 * \param[in] b first element index in the subtree
	 * \param[in] e one position past the last element index in the subtree
	 * \param[in] get_bbox a function that computes the bbox of an element
	 * \tparam GET_BBOX a function (or a functor) with the following arguments:
	 *  - mesh: a const reference to the mesh
	 *  - box: a reference where the computed bounding box of the element
	 *   will be stored
	 *  - element: the index of the element
	 */
	template<class GET_BBOX>
	void init_bboxes_recursive( const Mesh& M, std::vector<Box>& bboxes, index_t node_index, index_t b, index_t e, const GET_BBOX& get_bbox )
	{
		assert( node_index < bboxes.size() );
		assert( b != e );
		if( b + 1 == e )
		{
			get_bbox( M, bboxes[ node_index ], b );
			return;
		}
		index_t m = b + ( e - b ) / 2;
		index_t childl = 2 * node_index;
		index_t childr = 2 * node_index + 1;
		assert( childl < bboxes.size() );
		assert( childr < bboxes.size() );
		init_bboxes_recursive( M, bboxes, childl, b, m, get_bbox );
		init_bboxes_recursive( M, bboxes, childr, m, e, get_bbox );
		assert( childl < bboxes.size() );
		assert( childr < bboxes.size() );
		bbox_union( bboxes[ node_index ], bboxes[ childl ], bboxes[ childr ] );
	}

	inline void getFacetBox( const Mesh& mesh, Box32& box, uint32_t triangleIndex )
	{
		const vec3 *p1, *p2, *p3;
		mesh.getTriangleVertices( triangleIndex, &p1, &p2, &p3 );

		using namespace AvxMath;
		__m256d a = loadDouble3( &p1->x );
		__m256d b = loadDouble3( &p2->x );
		__m256d c = loadDouble3( &p3->x );

		box.computeTriangle( a, b, c );
	}

	static void initBoxesRecursive( const Mesh& mesh, std::vector<Box32>& bboxes, uint32_t node_index, index_t b, index_t e )
	{
		assert( node_index < bboxes.size() );
		assert( b != e );

		if( b + 1 == e )
		{
			getFacetBox( mesh, bboxes[ node_index ], b );
			return;
		}

		uint32_t m = b + ( e - b ) / 2;
		uint32_t childl = 2 * node_index;
		uint32_t childr = 2 * node_index + 1;
		assert( childl < bboxes.size() );
		assert( childr < bboxes.size() );
		initBoxesRecursive( mesh, bboxes, childl, b, m );
		initBoxesRecursive( mesh, bboxes, childr, m, e );
		bboxes[ node_index ].computeUnion( bboxes[ childl ], bboxes[ childr ] );
	}

	/**
	 * \brief Computes the squared distance between a point and a Box
	 *  with negative sign if the point is inside the Box.
	 * \param[in] p the point
	 * \param[in] B the box
	 * \return the signed squared distance between \p p and \p B
	 */

	static double point_box_signed_squared_distance_orig( const vec3& p, const Box& B )
	{
		bool inside = true;
		double result = 0.0;
		for( coord_index_t c = 0; c < 3; c++ )
		{
			if( p[ c ] < B.xyz_min[ c ] )
			{
				inside = false;
				result += geo_sqr( p[ c ] - B.xyz_min[ c ] );
			}
			else if( p[ c ] > B.xyz_max[ c ] )
			{
				inside = false;
				result += geo_sqr( p[ c ] - B.xyz_max[ c ] );
			}
		}
		if( inside )
		{
			result = -inner_point_box_squared_distance( p, B );
		}
		return result;
	}

#ifndef __AVX__
#error This code requires at least AVX1
#endif
	static __m128d pointBoxSignedSquaredDistance2( const __m256d pos, const Box& B )
	{
		using namespace AvxMath;
		// Load the box into 2 vectors
		const __m256d boxMin = _mm256_loadu_pd( &B.xyz_min[ 0 ] );
		const __m256d boxMax = _mm256_loadu_pd( &B.xyz_max[ 0 ] );

		// When inside, both numbers are positive
		const __m256d dmin = _mm256_sub_pd( pos, boxMin );
		const __m256d dmax = _mm256_sub_pd( boxMax, pos );
		// When inside, distance to the box
		// When outside, one of the vectors was negative another one positive, min will return the negative one, which is the distance we're after
		__m256d dist = _mm256_min_pd( dmin, dmax );

		// Compute vector mask of negative lanes, these are the lanes which were outside the box
		const __m256d zero = _mm256_setzero_pd();
		const __m256d outsideMaskVector = _mm256_cmp_pd( dist, zero, _CMP_LT_OQ );

		__m128d resIn, resOut;
		{
			// Compute negative distance for the inside case
			__m128d xy = low2( dist );
			__m128d z = high2( dist );
			xy = _mm_min_sd( xy, _mm_unpackhi_pd( xy, xy ) );
			xy = _mm_min_sd( xy, z );
			// Compute square
			xy = _mm_mul_sd( xy, xy );
			// Duplicate the low lane
			xy = _mm_movedup_pd( xy );
			// Negate the vector
			resIn = _mm_sub_pd( _mm256_castpd256_pd128( zero ), xy );
		}

		{
			// Compute positive distance for the outside case
			__m256d distOut = _mm256_and_pd( dist, outsideMaskVector );
			resOut = vector3Dot2( distOut, distOut );
		}

		// Compute bitwise OR of first 3 lanes in the outsideMaskVector, broadcasted to 2 lanes of the SSE vector
		__m128d maskHigh = _mm256_extractf128_pd( outsideMaskVector, 1 );
		__m128d maskLow = _mm256_castpd256_pd128( outsideMaskVector );
		maskLow = _mm_or_pd( maskLow, _mm_permute_pd( maskLow, _MM_SHUFFLE2( 0, 1 ) ) );
		maskLow = _mm_or_pd( maskLow, _mm_movedup_pd( maskHigh ) );

		// Select one of the output vectors
		return _mm_blendv_pd( resIn, resOut, maskLow );
	}

	static double pointBoxSignedSquaredDistance( const __m256d pos, const Box& B )
	{
		using namespace AvxMath;
		// Load into 3 vectors
		const __m256d boxMin = _mm256_loadu_pd( &B.xyz_min[ 0 ] );
		const __m256d boxMax = loadDouble3( &B.xyz_max[ 0 ] );

		// When inside, both numbers are positive
		const __m256d dmin = _mm256_sub_pd( pos, boxMin );
		const __m256d dmax = _mm256_sub_pd( boxMax, pos );
		// When inside, distance to the box; when outside, one of the vectors negative another positive, min will return the negative one
		__m256d dist = _mm256_min_pd( dmin, dmax );
		if( 0 == ( _mm256_movemask_pd( dist ) & 0b111 ) )
		{
			// The XYZ lanes of the dist are all non-negative, which means the point is inside the box
			// Compute horizontal minimum of the dist.xyz vector
			__m128d xy = low2( dist );
			__m128d z = high2( dist );
			xy = _mm_min_sd( xy, _mm_unpackhi_pd( xy, xy ) );
			xy = _mm_min_sd( xy, z );
			// Compute square
			xy = _mm_mul_sd( xy, xy );
			// Negate
			xy = _mm_sub_pd( _mm_setzero_pd(), xy );
			return _mm_cvtsd_f64( xy );
		}
		else
		{
			// Compute vector mask of negative lanes, these are the lanes which were outside the box
			__m256d outsideMaskVector = _mm256_cmp_pd( dist, _mm256_setzero_pd(), _CMP_LT_OQ );
			// Zero out the lanes for coordinates which were inside the box
			dist = _mm256_and_pd( dist, outsideMaskVector );
			// The square distance to the box is the length of that vector
			return vector3DotScalar( dist, dist );
		}
	}

	static double pointBoxSignedSquaredDistance( const vec3& p, const Box& B )
	{
		__m256d pos = AvxMath::loadDouble3( &p.x );
		return pointBoxSignedSquaredDistance( pos, B );
	}
}  // namespace

/****************************************************************************/

namespace floatTetWild
{
	MeshFacetsAABBWithEps::MeshFacetsAABBWithEps( const Mesh& M, std::vector<FacetRecursionStack>& stacks )
		: mesh_( M )
		, recursionStacks( stacks )
	{
		const uint32_t countTriangles = mesh_.countTriangles();
		const size_t boxesCount = max_node_index( 1, 0, countTriangles ) + 1;  //< this is because size == max_index + 1
		// Reserving one extra box so we can use full-vector unaligned loads to load coordinates from boxes, and not crash with access violations

		if constexpr( useFp32Boxes )
		{
			boxesFloat.reserve( boxesCount + 1 );
			boxesFloat.resize( boxesCount );
			initBoxesRecursive( M, boxesFloat, 1, 0, countTriangles );
		}

		if constexpr( !useFp32Boxes || dbgCompareVersions )
		{
			bboxes_.reserve( boxesCount + 1 );
			bboxes_.resize( boxesCount );
			init_bboxes_recursive( mesh_, bboxes_, 1, 0, countTriangles, get_facet_bbox );
		}
	}

	void MeshFacetsAABBWithEps::nearestFacetRecursive(
	  const vec3& p, index_t& nearest_f, vec3& nearest_point, double& sq_dist, index_t n, index_t b, index_t e ) const
	{
		assert( e > b );

		// If node is a leaf: compute point-facet distance and replace current if nearer
		const __m256d pos = AvxMath::loadDouble3( &p.x );
		if( b + 1 == e )
		{
			vec3 cur_nearest_point;
			double cur_sq_dist;
			get_point_facet_nearest_point( mesh_, pos, b, cur_nearest_point, cur_sq_dist );
			if( cur_sq_dist < sq_dist )
			{
				nearest_f = b;
				nearest_point = cur_nearest_point;
				sq_dist = cur_sq_dist;
			}
			return;
		}
		index_t m = b + ( e - b ) / 2;
		index_t childl = 2 * n;
		index_t childr = 2 * n + 1;

		double dl = pointBoxSignedSquaredDistance( pos, bboxes_[ childl ] );
		double dr = pointBoxSignedSquaredDistance( pos, bboxes_[ childr ] );

		// Traverse the "nearest" child first, so that it has more chances
		// to prune the traversal of the other child.
		if( dl < dr )
		{
			if( dl < sq_dist )
			{
				nearestFacetRecursive( p, nearest_f, nearest_point, sq_dist, childl, b, m );
			}
			if( dr < sq_dist )
			{
				nearestFacetRecursive( p, nearest_f, nearest_point, sq_dist, childr, m, e );
			}
		}
		else
		{
			if( dr < sq_dist )
			{
				nearestFacetRecursive( p, nearest_f, nearest_point, sq_dist, childr, m, e );
			}
			if( dl < sq_dist )
			{
				nearestFacetRecursive( p, nearest_f, nearest_point, sq_dist, childl, b, m );
			}
		}
	}

	void MeshFacetsAABBWithEps::facetInEnvelopeRecursive(
	  __m256d p, double sq_epsilon, index_t& nearest_f, vec3& nearest_point, double& sq_dist, __m128i nbe ) const
	{
		const uint32_t n = (uint32_t)_mm_cvtsi128_si32( nbe );
		const uint32_t b = (uint32_t)_mm_extract_epi32( nbe, 1 );
		const uint32_t e = (uint32_t)_mm_extract_epi32( nbe, 2 );
		assert( e > b );

		if( sq_dist <= sq_epsilon )
			return;

		// If node is a leaf: compute point-facet distance
		// and replace current if nearer
		if( b + 1 == e )
		{
			vec3 cur_nearest_point;
			double cur_sq_dist;
			get_point_facet_nearest_point( mesh_, p, b, cur_nearest_point, cur_sq_dist );
			if( cur_sq_dist < sq_dist )
			{
				nearest_f = b;
				nearest_point = cur_nearest_point;
				sq_dist = cur_sq_dist;
			}
			return;
		}
		index_t m = b + ( e - b ) / 2;
		index_t childl = 2 * n;
		index_t childr = 2 * n + 1;

		// The original code suffers from the unpredictable "if( dl < dr )" branch
		// To workaround, we using vector blends as conditional moves to figure out which way to go first.
		const __m128d dl = pointBoxSignedSquaredDistance2( p, bboxes_[ childl ] );
		const __m128d dr = pointBoxSignedSquaredDistance2( p, bboxes_[ childr ] );

		// Compare for dl < dr
		// Because pointBoxSignedSquaredDistance2 returns vectors with both lanes equal, the result is either zero or a vector of UINT_MAX
		const __m128i lt = _mm_castpd_si128( _mm_cmplt_pd( dl, dr ) );

		// Create left/right index vectors
		const __m128i recLeft = _mm_setr_epi32( (int)childl, (int)b, (int)m, 0 );
		const __m128i recRight = _mm_setr_epi32( (int)childr, (int)m, (int)e, 0 );

		double d = _mm_cvtsd_f64( _mm_min_sd( dr, dl ) );
		// The remaining 2 branches are very predictable, they are both almost always taken
		if( d < sq_dist && d <= sq_epsilon )
		{
			const __m128i rec = _mm_blendv_epi8( recRight, recLeft, lt );  // ( dl < dr ) ? recLeft : recRight
			facetInEnvelopeRecursive( p, sq_epsilon, nearest_f, nearest_point, sq_dist, rec );
		}

		d = _mm_cvtsd_f64( _mm_max_sd( dr, dl ) );
		if( d < sq_dist && d <= sq_epsilon )
		{
			const __m128i rec = _mm_blendv_epi8( recLeft, recRight, lt );  // ( dl < dr ) ? recRight : recLeft
			facetInEnvelopeRecursive( p, sq_epsilon, nearest_f, nearest_point, sq_dist, rec );
		}
	}

	inline __m128i makeUInt3( uint32_t x, uint32_t y, uint32_t z )
	{
		__m128i v = _mm_cvtsi32_si128( (int)x );
		v = _mm_insert_epi32( v, (int)y, 1 );
		v = _mm_insert_epi32( v, (int)z, 2 );
		return v;
	}

	inline void unpackUInt3( __m128i vec, uint32_t& x, uint32_t& y, uint32_t& z )
	{
		x = (uint32_t)_mm_cvtsi128_si32( vec );
		y = (uint32_t)_mm_extract_epi32( vec, 1 );
		z = (uint32_t)_mm_extract_epi32( vec, 2 );
	}

	void MeshFacetsAABBWithEps::nearestFacetStack( const vec3& point, index_t& nearestFacet, vec3& nearestPoint, double& sqDistResult ) const
	{
		std::vector<FacetRecursionFrame>& stack = recursionStacks[ omp_get_thread_num() ].stack;
		const __m256d p = AvxMath::loadDouble3( &point.x );

		// Copy that value from memory to a register, saves quite a few loads/stores in the loop below
		double sqDist = sqDistResult;

		// Setup the initial state
		uint32_t n = 1;
		uint32_t b = 0;
		uint32_t e = (uint32_t)mesh_.countTriangles();
		double d = 0.0;

#define POP_FROM_THE_STACK()                     \
	if( stack.empty() )                          \
		break;                                   \
	const FacetRecursionFrame& f = stack.back(); \
	n = f.n;                                     \
	b = f.b;                                     \
	e = f.e;                                     \
	d = f.d;                                     \
	stack.pop_back()

		// Run the "recursion" using an std::vector instead of the stack
		while( true )
		{
			assert( e > b );

			if( d >= sqDist )
			{
				// The original version would have skipped this frame with "if( d < sq_dist )"
				POP_FROM_THE_STACK();
				continue;
			}

			if( b + 1 == e )
			{
				// If node is a leaf: compute point-facet distance and replace current if nearer
				vec3 cur_nearest_point;
				double cur_sq_dist;
				get_point_facet_nearest_point( mesh_, p, b, cur_nearest_point, cur_sq_dist );
				if( cur_sq_dist < sqDist )
				{
					nearestFacet = b;
					nearestPoint = cur_nearest_point;
					sqDist = cur_sq_dist;
				}
				POP_FROM_THE_STACK();
				continue;
			}

			const uint32_t m = b + ( e - b ) / 2;
			const uint32_t childl = 2 * n;
			const uint32_t childr = 2 * n + 1;

			// The original code suffers from the unpredictable "if( dl < dr )" branch
			// To workaround, we using vector blends as conditional moves to figure out which way to go first.
			const __m128d dl = pointBoxSignedSquaredDistance2( p, bboxes_[ childl ] );
			const __m128d dr = pointBoxSignedSquaredDistance2( p, bboxes_[ childr ] );

			// Compare for dl < dr
			// Because pointBoxSignedSquaredDistance2 returns vectors with both lanes equal, the result is either zero or a vector of UINT_MAX
			const __m128i lt = _mm_castpd_si128( _mm_cmplt_pd( dl, dr ) );

			// Create left/right index vectors
			const __m128i recLeft = makeUInt3( childl, b, m );
			const __m128i recRight = makeUInt3( childr, m, e );

			// Push the SECOND recursive call of the original version to the stack
			__m128i rec = _mm_blendv_epi8( recLeft, recRight, lt );
			FacetRecursionFrame& newFrame = stack.emplace_back();
			newFrame.storeIndices( rec );
			newFrame.d = _mm_cvtsd_f64( _mm_max_sd( dr, dl ) );

			// Replace the local variables with the FIRST recursive call of the original version
			rec = _mm_blendv_epi8( recRight, recLeft, lt );
			unpackUInt3( rec, n, b, e );
			d = _mm_cvtsd_f64( _mm_min_sd( dr, dl ) );
		}
#undef POP_FROM_THE_STACK

		assert( stack.empty() );
		// Store the result back to memory
		sqDistResult = sqDist;
	}

	void MeshFacetsAABBWithEps::nearestFacetStack32( const vec3& point, index_t& nearestFacet, vec3& nearestPoint, double& sqDistResult ) const
	{
		std::vector<FacetRecursionFrame32>& stack = recursionStacks[ omp_get_thread_num() ].stack32;
		const __m256d p = AvxMath::loadDouble3( &point.x );
		const __m128 p32 = _mm256_cvtpd_ps( p );

		// Copy that value from memory to a register, saves quite a few loads/stores in the loop below
		double sqDist = sqDistResult;

		// Setup the initial state
		uint32_t n = 1;
		uint32_t b = 0;
		uint32_t e = (uint32_t)mesh_.countTriangles();
		float d = 0.0f;

#define POP_FROM_THE_STACK()                       \
	if( stack.empty() )                            \
		break;                                     \
	const FacetRecursionFrame32& f = stack.back(); \
	n = f.n;                                       \
	b = f.b;                                       \
	e = f.e;                                       \
	d = f.d;                                       \
	stack.pop_back()

		// Run the "recursion" using an std::vector instead of the stack
		while( true )
		{
			assert( e > b );

			if( d >= sqDist )
			{
				// The original version would have skipped this frame with "if( d < sq_dist )"
				POP_FROM_THE_STACK();
				continue;
			}

			if( b + 1 == e )
			{
				// If node is a leaf: compute point-facet distance and replace current if nearer
				vec3 cur_nearest_point;
				double cur_sq_dist;
				get_point_facet_nearest_point( mesh_, p, b, cur_nearest_point, cur_sq_dist );
				if( cur_sq_dist < sqDist )
				{
					nearestFacet = b;
					nearestPoint = cur_nearest_point;
					sqDist = cur_sq_dist;
				}
				POP_FROM_THE_STACK();
				continue;
			}

			const uint32_t m = b + ( e - b ) / 2;
			const uint32_t childl = 2 * n;
			const uint32_t childr = 2 * n + 1;

			// The original code suffers from the unpredictable "if( dl < dr )" branch
			// To workaround, we using vector blends as conditional moves to figure out which way to go first.
			const __m128 dl = boxesFloat[ childl ].pointBoxSignedSquaredDistance( p32 );
			const __m128 dr = boxesFloat[ childr ].pointBoxSignedSquaredDistance( p32 );

			// Compare for dl < dr
			// Because pointBoxSignedSquaredDistance2 returns vectors with both lanes equal, the result is either zero or a vector of UINT_MAX
			const __m128i lt = _mm_castps_si128( _mm_cmplt_ps( dl, dr ) );

			// Create left/right index vectors
			const __m128i recLeft = makeUInt3( childl, b, m );
			const __m128i recRight = makeUInt3( childr, m, e );

			// Push the SECOND recursive call of the original version to the stack
			__m128i rec = _mm_blendv_epi8( recLeft, recRight, lt );
			FacetRecursionFrame32& newFrame = stack.emplace_back();
			newFrame.store( rec, _mm_max_ps( dr, dl ) );

			// Replace the local variables with the FIRST recursive call of the original version
			rec = _mm_blendv_epi8( recRight, recLeft, lt );
			unpackUInt3( rec, n, b, e );
			d = _mm_cvtss_f32( _mm_min_ss( dr, dl ) );
		}
#undef POP_FROM_THE_STACK

		assert( stack.empty() );
		// Store the result back to memory
		sqDistResult = sqDist;
	}

	void MeshFacetsAABBWithEps::facetInEnvelopeStack(
	  __m256d p, double sqEpsilon, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDistResult ) const
	{
		std::vector<FacetRecursionFrame>& stack = recursionStacks[ omp_get_thread_num() ].stack;

		// Copy that value from memory to a register, saves quite a few loads/stores in the loop below
		double sqDist = sqDistResult;

		// Setup the initial state
		uint32_t n = 1;
		uint32_t b = 0;
		uint32_t e = (uint32_t)mesh_.countTriangles();
		double d = 0.0;

#define POP_FROM_THE_STACK()                     \
	if( stack.empty() )                          \
		break;                                   \
	const FacetRecursionFrame& f = stack.back(); \
	n = f.n;                                     \
	b = f.b;                                     \
	e = f.e;                                     \
	d = f.d;                                     \
	stack.pop_back()

		// Run the "recursion" using an std::vector instead of the stack
		while( true )
		{
			assert( e > b );

			if( d >= sqDist || d > sqEpsilon )
			{
				// The original version would have skipped this frame with "if( d < sq_dist && d <= sq_epsilon )"
				POP_FROM_THE_STACK();
				continue;
			}

			if( b + 1 == e )
			{
				// If node is a leaf: compute point-facet distance and replace current if nearer
				vec3 cur_nearest_point;
				double cur_sq_dist;
				get_point_facet_nearest_point( mesh_, p, b, cur_nearest_point, cur_sq_dist );
				if( cur_sq_dist < sqDist )
				{
					nearestFacet = b;
					nearestPoint = cur_nearest_point;
					sqDist = cur_sq_dist;
					if( cur_sq_dist <= sqEpsilon )
					{
						// This alone saves non-trivial amount of overhead compared to recursive version
						// Clearing the complete std::vector only takes a few instructions
						stack.clear();
						break;
					}
				}
				POP_FROM_THE_STACK();
				continue;
			}

			const uint32_t m = b + ( e - b ) / 2;
			const uint32_t childl = 2 * n;
			const uint32_t childr = 2 * n + 1;

			// The original code suffers from the unpredictable "if( dl < dr )" branch
			// To workaround, we using vector blends as conditional moves to figure out which way to go first.
			const __m128d dl = pointBoxSignedSquaredDistance2( p, bboxes_[ childl ] );
			const __m128d dr = pointBoxSignedSquaredDistance2( p, bboxes_[ childr ] );

			// Compare for dl < dr
			// Because pointBoxSignedSquaredDistance2 returns vectors with both lanes equal, the result is either zero or a vector of UINT_MAX
			const __m128i lt = _mm_castpd_si128( _mm_cmplt_pd( dl, dr ) );

			// Create left/right index vectors
			const __m128i recLeft = makeUInt3( childl, b, m );
			const __m128i recRight = makeUInt3( childr, m, e );

			// Push the SECOND recursive call of the original version to the stack
			__m128i rec = _mm_blendv_epi8( recLeft, recRight, lt );
			FacetRecursionFrame& newFrame = stack.emplace_back();
			newFrame.storeIndices( rec );
			newFrame.d = _mm_cvtsd_f64( _mm_max_sd( dr, dl ) );

			// Replace the local variables with the FIRST recursive call of the original version
			rec = _mm_blendv_epi8( recRight, recLeft, lt );
			unpackUInt3( rec, n, b, e );
			d = _mm_cvtsd_f64( _mm_min_sd( dr, dl ) );
		}
#undef POP_FROM_THE_STACK

		assert( stack.empty() );
		// Store the result back to memory
		sqDistResult = sqDist;
	}

	void MeshFacetsAABBWithEps::facetInEnvelopeStack32(
	  __m256d p, double sqEpsilon, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDistResult ) const
	{
		std::vector<FacetRecursionFrame32>& stack = recursionStacks[ omp_get_thread_num() ].stack32;

		// Copy that value from memory to a register, saves quite a few loads/stores in the loop below
		double sqDist = sqDistResult;

		// Setup the initial state
		uint32_t n = 1;
		uint32_t b = 0;
		uint32_t e = (uint32_t)mesh_.countTriangles();
		float d = 0.0f;
		const __m128 p32 = _mm256_cvtpd_ps( p );

#define POP_FROM_THE_STACK()      \
	if( stack.empty() )           \
		break;                    \
	const auto& f = stack.back(); \
	n = f.n;                      \
	b = f.b;                      \
	e = f.e;                      \
	d = f.d;                      \
	stack.pop_back()

		// Run the "recursion" using an std::vector instead of the stack
		while( true )
		{
			assert( e > b );

			if( d >= sqDist || d > sqEpsilon )
			{
				// The original version would have skipped this frame with "if( d < sq_dist && d <= sq_epsilon )"
				POP_FROM_THE_STACK();
				continue;
			}

			if( b + 1 == e )
			{
				// If node is a leaf: compute point-facet distance and replace current if nearer
				vec3 cur_nearest_point;
				double cur_sq_dist;
				get_point_facet_nearest_point( mesh_, p, b, cur_nearest_point, cur_sq_dist );
				if( cur_sq_dist < sqDist )
				{
					nearestFacet = b;
					nearestPoint = cur_nearest_point;
					sqDist = cur_sq_dist;
					if( cur_sq_dist <= sqEpsilon )
					{
						// This alone saves non-trivial amount of overhead compared to recursive version
						// Clearing the complete std::vector only takes a few instructions
						stack.clear();
						break;
					}
				}
				POP_FROM_THE_STACK();
				continue;
			}

			const uint32_t m = b + ( e - b ) / 2;
			const uint32_t childl = 2 * n;
			const uint32_t childr = 2 * n + 1;

			// The original code suffers from the unpredictable "if( dl < dr )" branch
			// To workaround, we using vector blends as conditional moves to figure out which way to go first.
			const __m128 dl = boxesFloat[ childl ].pointBoxSignedSquaredDistance( p32 );
			const __m128 dr = boxesFloat[ childr ].pointBoxSignedSquaredDistance( p32 );

			// Compare for dl < dr
			// Because pointBoxSignedSquaredDistance2 returns vectors with both lanes equal, the result is either zero or a vector of UINT_MAX
			const __m128i lt = _mm_castps_si128( _mm_cmplt_ps( dl, dr ) );

			// Create left/right index vectors
			const __m128i recLeft = makeUInt3( childl, b, m );
			const __m128i recRight = makeUInt3( childr, m, e );

			// Push the SECOND recursive call of the original version to the stack
			__m128i rec = _mm_blendv_epi8( recLeft, recRight, lt );
			auto& newFrame = stack.emplace_back();
			newFrame.store( rec, _mm_max_ps( dr, dl ) );

			// Replace the local variables with the FIRST recursive call of the original version
			rec = _mm_blendv_epi8( recRight, recLeft, lt );
			unpackUInt3( rec, n, b, e );
			d = _mm_cvtss_f32( _mm_min_ss( dr, dl ) );
		}
#undef POP_FROM_THE_STACK

		assert( stack.empty() );
		// Store the result back to memory
		sqDistResult = sqDist;
	}

	void MeshFacetsAABBWithEps::facetInEnvelopeCompare(
	  __m256d p, double sqEpsilon, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDist, __m128i nbe ) const
	{
		const uint32_t nfInput = nearestFacet;
		const vec3 npInput = nearestPoint;
		const double sqDistInput = sqDist;

		uint32_t nf1 = nfInput;
		vec3 np1 = npInput;
		double sqd1 = sqDistInput;
		facetInEnvelopeRecursive( p, sqEpsilon, nf1, np1, sqd1, nbe );

		nearestFacet = nfInput;
		nearestPoint = npInput;
		sqDist = sqDistInput;
		facetInEnvelopeStack( p, sqEpsilon, nearestFacet, nearestPoint, sqDist );

		if( sqDist <= sqd1 )
		{
			// Despite I tried to replicate the original recursive version, the new version sometimes finds closer points than the old one
			// I have no idea why, let's hope that's a good thing
			return;
		}

		__debugbreak();
	}

	void MeshFacetsAABBWithEps::nearestFacetCompare( const GEO2::vec3& p, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDist ) const
	{
		const uint32_t nfInput = nearestFacet;
		const vec3 npInput = nearestPoint;
		const double sqDistInput = sqDist;

		uint32_t nf1 = nfInput;
		vec3 np1 = npInput;
		double sqd1 = sqDistInput;
		nearestFacetRecursive( p, nf1, np1, sqd1, 1, 0, mesh_.countTriangles() );

		nearestFacet = nfInput;
		nearestPoint = npInput;
		sqDist = sqDistInput;
		nearestFacetStack( p, nearestFacet, nearestPoint, sqDist );

		if( sqDist == sqd1 )
			return;
		__debugbreak();
	}
}  // namespace floatTetWild
