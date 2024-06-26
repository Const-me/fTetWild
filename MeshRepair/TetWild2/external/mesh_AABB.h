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

#pragma once
/**
 * \file mesh_AABB.h
 * \brief Axis Aligned Bounding Box trees for accelerating
 *  geometric queries that operate on a Mesh.
 */

#include <string>
#include "../Utils/Geogram2.h"
#include "../Utils/Mesh/TriangleMesh.h"
#include "../src/FacetRecursionStack.h"
#include "../Utils/Box32.h"

namespace floatTetWild
{
	/**
	 * \brief Axis Aligned Bounding Box tree of mesh facets.
	 * \details Used to quickly compute facet intersection and
	 *  to locate the nearest facet from 3d query points.
	 */
	class MeshFacetsAABBWithEps
	{
	  public:
		/**
		 * \brief Creates the Axis Aligned Bounding Boxes tree.
		 * \param[in] M the input mesh. It can be modified,
		 *  and will be triangulated (if
		 *  not already a triangular mesh). The facets are
		 *  re-ordered (using Morton's order, see mesh_reorder()).
		 * \param[in] reorder if not set, Morton re-ordering is
		 *  skipped (but it means that mesh_reorder() was previously
		 *  called else the algorithm will be pretty unefficient).
		 * \pre M.facets.are_simplices()
		 */
		MeshFacetsAABBWithEps( const GEO2::Mesh& M, FacetRecursionStacks& stacks );

		/**
		 * \brief Finds the nearest facet from an arbitrary 3d query point.
		 * \param[in] p query point
		 * \param[out] nearest_point nearest point on the surface
		 * \param[out] sq_dist squared distance between p and the surface.
		 * \return the index of the facet nearest to point p.
		 */
		GEO2::index_t nearest_facet( const GEO2::vec3& p, GEO2::vec3& nearest_point, double& sq_dist ) const
		{
			GEO2::index_t nearest_facet = GEO2::NO_FACET;
			sq_dist = DBL_MAX;
			nearestFacetStack32( p, nearest_facet, nearest_point, sq_dist );
			return nearest_facet;
		}

		/**
		 * \brief Computes the nearest point and nearest facet from
		 * a query point, using user-specified hint.
		 *
		 * \details The hint is specified as reasonable initial values of
		 * (nearest_facet, nearest_point, sq_dist). If multiple queries
		 * are done on a set of points that has spatial locality,
		 * the hint can be the result of the previous call.
		 *
		 * \param[in] p query point
		 * \param[in,out] nearest_facet the nearest facet so far,
		 *   or NO_FACET if not known yet
		 * \param[in,out] nearest_point a point in nearest_facet
		 * \param[in,out] sq_dist squared distance between p and
		 *    nearest_point
		 * \note On entry, \p sq_dist needs to be equal to the squared
		 *   distance between \p p and \p nearest_point (it is easy to
		 *   forget to update it when calling it within a loop).
		 */
		void nearest_facet_with_hint( const GEO2::vec3& p, GEO2::index_t& nearest_facet, GEO2::vec3& nearest_point, double& sq_dist ) const
		{
			if( nearest_facet == GEO2::NO_FACET )
				sq_dist = DBL_MAX;
			nearestFacetStack32( p, nearest_facet, nearest_point, sq_dist );
		}

		/*
		 * Finds the nearest facet on the surface, but stops early if a
		 * point within a given distance is found.
		 */
		GEO2::index_t facet_in_envelope( const GEO2::vec3& p, double sq_epsilon, GEO2::vec3& nearest_point, double& sq_dist ) const
		{
			GEO2::index_t nearest_facet = GEO2::NO_FACET;
			sq_dist = DBL_MAX;
			__m256d pt = AvxMath::loadDouble3( &p.x );
			facetInEnvelopeStack32( pt, sq_epsilon, nearest_facet, nearest_point, sq_dist );
			return nearest_facet;
		}

		/*
		 * Same as before, but stops as soon as a point on the surface in
		 * within a given distance bound from the triangle mesh.
		 */
		void facet_in_envelope_with_hint(
		  const GEO2::vec3& p, double sq_epsilon, GEO2::index_t& nearest_facet, GEO2::vec3& nearest_point, double& sq_dist ) const
		{
			if( nearest_facet == GEO2::NO_FACET )
				sq_dist = DBL_MAX;
			__m256d pt = AvxMath::loadDouble3( &p.x );
			facetInEnvelopeStack32( pt, sq_epsilon, nearest_facet, nearest_point, sq_dist );
		}

		// Find root node of the tree which contains the specified box
		__m128i __vectorcall getBoxRoot( __m256d boxMin64, __m256d boxMax64 ) const;

		bool __vectorcall isOutOfEnvelope( __m256d pos, __m128d eps21, __m128i searchRoot, uint32_t& prevFace ) const;
		bool __vectorcall isOutOfEnvelope( __m256d pos, __m128d eps21, uint32_t& prevFace ) const;

		// Collect IDs of the faces which might intersect the specified bounding box
		// The returned vector is stored in a thread-local structure for the calling thread
		const std::vector<uint32_t>& facesInTheBox( __m256d boxMin64, __m256d boxMax64 ) const;

		// Compute minimum squared distance between the point and the specified set of triangles
		// Return true if that distance exceeds the supplied scalar
		bool isOutOfEnvelope( __m256d pos, double eps2, const std::vector<uint32_t>& faces ) const;

	  private:

		// Same as above without recursion
		void nearestFacetStack32( const GEO2::vec3& point, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDistResult ) const;

		void facetInEnvelopeStack32_v1( __m256d p, double sqEpsilon, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDist ) const;
		void facetInEnvelopeStack32_v2( __m256d p, double sqEpsilon, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDist ) const;
		void facetInEnvelopeStack32( __m256d p, double sqEpsilon, GEO2::index_t& nearestFacet, GEO2::vec3& nearestPoint, double& sqDist ) const;

		std::vector<Box32> boxesFloat;
		const GEO2::Mesh& mesh_;
		FacetRecursionStacks& recursionStacks;
	};

	inline void get_point_facet_nearest_point( const GEO2::Mesh& M, __m256d p, GEO2::index_t f, GEO2::vec3& nearest_p, double& squared_dist )
	{
		using namespace GEO2;
		const vec3 *p1, *p2, *p3;
		M.getTriangleVertices( f, &p1, &p2, &p3 );
		squared_dist = point_triangle_squared_distance( p, *p1, *p2, *p3, &nearest_p );
	}

	inline double pointTriangleSquaredDistance( const GEO2::Mesh& M, __m256d p, GEO2::index_t f )
	{
		using namespace GEO2;
		const vec3 *p1, *p2, *p3;
		M.getTriangleVertices( f, &p1, &p2, &p3 );
		return point_triangle_squared_distance( p, *p1, *p2, *p3, nullptr );
	}
}  // namespace floatTetWild