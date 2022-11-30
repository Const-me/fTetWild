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
#include <Utils/AvxMath.h>
#include <Utils/miscUtils.h>
#include <omp.h>

namespace
{
	using namespace GEO2;
	using namespace floatTetWild;

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

	// Compute axis-aligned bounding box of a mesh facet, and downcast it to FP32 with proper rounding: towards -INF for minimum, and towards +INF for maximum
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

	// Compute the hierarchy of bounding boxes recursively
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

		boxesFloat.resize( boxesCount + 1 );
		boxesFloat.back().setZero();
		boxesFloat.pop_back();

		initBoxesRecursive( M, boxesFloat, 1, 0, countTriangles );
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

	void MeshFacetsAABBWithEps::nearestFacetStack32( const vec3& point, index_t& nearestFacet, vec3& nearestPoint, double& sqDistResult ) const
	{
		std::vector<FacetRecursionFrame32>& stack = recursionStacks[ omp_get_thread_num() ].stack32;
		const __m256d p = AvxMath::loadDouble3( &point.x );
		const __m128 p32_single = _mm256_cvtpd_ps( p );
		// Duplicate the downcasted position into low & high halves of an AVX vector, for pointBoxSignedSquaredDistanceX2 function
		const __m256 p32 = _mm256_setr_m128( p32_single, p32_single );

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

			const __m256 distX2 = Box32::pointBoxSignedSquaredDistanceX2( boxesFloat[ childl ], boxesFloat[ childr ], p32 );
			const __m128 dl = _mm256_castps256_ps128( distX2 );
			const __m128 dr = _mm256_extractf128_ps( distX2, 1 );

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

	void MeshFacetsAABBWithEps::facetInEnvelopeStack32_v1(
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
		const __m128 p32_single = _mm256_cvtpd_ps( p );
		// Duplicate the downcasted position into low & high halves of an AVX vector, for pointBoxSignedSquaredDistanceX2 function
		const __m256 p32 = _mm256_setr_m128( p32_single, p32_single );

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

			const __m256 distX2 = Box32::pointBoxSignedSquaredDistanceX2( boxesFloat[ childl ], boxesFloat[ childr ], p32 );
#if 0
			// Testing things, compare against older pointBoxSignedSquaredDistance implementation
			const __m128 dl_old = boxesFloat[ childl ].pointBoxSignedSquaredDistance( _mm256_castps256_ps128( p32 ) );
			const __m128 dr_old = boxesFloat[ childr ].pointBoxSignedSquaredDistance( _mm256_castps256_ps128( p32 ) );
			const __m256 distX2_old = _mm256_setr_m128( dl_old, dr_old );
			const __m256 neq = _mm256_cmp_ps( distX2, distX2_old, _CMP_NEQ_UQ );
			if( !_mm256_testz_ps( neq, neq ) )
				__debugbreak();
#endif
			const __m128 dl = _mm256_castps256_ps128( distX2 );
			const __m128 dr = _mm256_extractf128_ps( distX2, 1 );

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

	inline __m256 makeBoxVector( __m256d p, double sqEpsilon )
	{
		// std::sqrt() is a function call, this code is called too frequently, compute square root with 1 instruction instead
		__m128d e2 = _mm_set_sd( sqEpsilon );
		__m128d eps = _mm_sqrt_sd( e2, e2 );

		const __m256d ev = _mm256_set1_pd( _mm_cvtsd_f64( eps ) );
		const __m256d i = _mm256_sub_pd( p, ev );
		const __m256d ax = _mm256_add_pd( p, ev );
		return Box32::createBoxVector( i, ax );
	}

	void MeshFacetsAABBWithEps::facetInEnvelopeStack32_v2(
	  __m256d p, double sqEpsilon, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDistResult ) const
	{
		std::vector<FacetRecursionFrame32>& stack = recursionStacks[ omp_get_thread_num() ].stack32;
		const __m256 boxVec = makeBoxVector( p, sqEpsilon );

		// Copy that value from memory to a register, saves quite a few loads/stores in the loop below
		double sqDist = sqDistResult;

		// Setup the initial state
		uint32_t n = 1;
		uint32_t b = 0;
		uint32_t e = (uint32_t)mesh_.countTriangles();

#define POP_FROM_THE_STACK()      \
	if( stack.empty() )           \
		break;                    \
	const auto& f = stack.back(); \
	n = f.n;                      \
	b = f.b;                      \
	e = f.e;                      \
	stack.pop_back()

		// Run the "recursion" using an std::vector instead of the stack
		while( true )
		{
			assert( e > b );

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

			uint8_t bitmap = boxesFloat[ childl ].intersects( boxVec ) ? 1 : 0;
			bitmap |= boxesFloat[ childr ].intersects( boxVec ) ? 2 : 0;

			if( 0 == bitmap )
			{
				// None of the 2 children intersected with the query box
				POP_FROM_THE_STACK();
				continue;
			}

			if( 1 == bitmap )
			{
				// Only the left child has intersected, replace the state with [ childl, b, m ]
				n = childl;
				e = m;
				continue;
			}

			if( 2 == bitmap )
			{
				// Only the right child has intersected, replace the state with [ childr, m, e ]
				n = childr;
				b = m;
				continue;
			}

			// Push the right [ childr, m, e ] state to the stack..
			auto& newFrame = stack.emplace_back();
			newFrame.n = childr;
			newFrame.b = m;
			newFrame.e = e;
			// ..and replace the current state with the left one
			n = childl;
			e = m;
		}
#undef POP_FROM_THE_STACK

		assert( stack.empty() );
		// Store the result back to memory
		sqDistResult = sqDist;
	}

	void MeshFacetsAABBWithEps::facetInEnvelopeStack32(
	  __m256d p, double sqEpsilon, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDist ) const
	{
#if 1
		facetInEnvelopeStack32_v2( p, sqEpsilon, nearestFacet, nearestPoint, sqDist );
#else
		const uint32_t nfInput = nearestFacet;
		const vec3 npInput = nearestPoint;
		const double sqDistInput = sqDist;

		uint32_t nf1 = nfInput;
		vec3 np1 = npInput;
		double sqd1 = sqDistInput;
		facetInEnvelopeStack32_v1( p, sqEpsilon, nf1, np1, sqd1 );

		nearestFacet = nfInput;
		nearestPoint = npInput;
		sqDist = sqDistInput;
		facetInEnvelopeStack32_v2( p, sqEpsilon, nearestFacet, nearestPoint, sqDist );

		const bool found1 = sqd1 <= sqEpsilon;
		const bool found2 = sqDist <= sqEpsilon;
		if( found1 == found2 )
			return;
		// clang-format off
		printf( "AABB query: eps2 = %g, squared distances [ %g, %g ], found = [ %s, %s ], facets = [ %i, %i ]\n",
			sqEpsilon, sqd1, sqDist, cstr( found1 ), cstr( found2 ), nf1, nearestFacet );
		// clang-format on
#endif
	}

	const std::vector<uint32_t>& MeshFacetsAABBWithEps::facesInTheBox( __m256d boxMin64, __m256d boxMax64 ) const
	{
		auto& frs = recursionStacks[ omp_get_thread_num() ];
		std::vector<FacetRecursionFrame32>& stack = frs.stack32;
		std::vector<uint32_t>& faces = frs.faces;
		faces.clear();

		const __m256 boxVec = Box32::createBoxVector( boxMin64, boxMax64 );
		// Setup the initial state
		uint32_t n = 1;
		uint32_t b = 0;
		uint32_t e = (uint32_t)mesh_.countTriangles();

#define POP_FROM_THE_STACK()      \
	if( stack.empty() )           \
		break;                    \
	const auto& f = stack.back(); \
	n = f.n;                      \
	b = f.b;                      \
	e = f.e;                      \
	stack.pop_back()

		// Run the "recursion" using an std::vector instead of the stack
		while( true )
		{
			assert( e > b );

			if( b + 1 == e )
			{
				// Node is a leaf, add the triangle to the output collection.
				// TODO: maybe a better filtering here, like plane versus aligned box, or plane versus sphere
				faces.push_back( b );
				POP_FROM_THE_STACK();
				continue;
			}

			const uint32_t m = b + ( e - b ) / 2;
			const uint32_t childl = 2 * n;
			const uint32_t childr = 2 * n + 1;

			uint8_t bitmap = boxesFloat[ childl ].intersects( boxVec ) ? 1 : 0;
			bitmap |= boxesFloat[ childr ].intersects( boxVec ) ? 2 : 0;

			if( 0 == bitmap )
			{
				// None of the 2 children intersected with the query box
				POP_FROM_THE_STACK();
				continue;
			}

			if( 1 == bitmap )
			{
				// Only the left child has intersected, replace the state with [ childl, b, m ]
				n = childl;
				e = m;
				continue;
			}

			if( 2 == bitmap )
			{
				// Only the right child has intersected, replace the state with [ childr, m, e ]
				n = childr;
				b = m;
				continue;
			}

			// They both intersected. Push the right [ childr, m, e ] state to the stack, and replace the current state with the left one
			auto& newFrame = stack.emplace_back();
			newFrame.n = childr;
			newFrame.b = m;
			newFrame.e = e;

			n = childl;
			e = m;
		}
#undef POP_FROM_THE_STACK

		return faces;
	}

	bool MeshFacetsAABBWithEps::isOutOfEnvelope( __m256d pos, double eps2, const std::vector<uint32_t>& faces ) const
	{
		for( uint32_t f : faces )
		{
			const vec3 *p1, *p2, *p3;
			mesh_.getTriangleVertices( f, &p1, &p2, &p3 );

			const double sqd = point_triangle_squared_distance( pos, *p1, *p2, *p3, nullptr );
			if( sqd < eps2 )
				return false;
		}
		return true;
	}
}  // namespace floatTetWild