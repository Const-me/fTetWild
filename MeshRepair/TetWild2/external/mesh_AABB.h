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
		MeshFacetsAABBWithEps( const GEO2::Mesh& M );

		/**
		 * \brief Finds the nearest facet from an arbitrary 3d query point.
		 * \param[in] p query point
		 * \param[out] nearest_point nearest point on the surface
		 * \param[out] sq_dist squared distance between p and the surface.
		 * \return the index of the facet nearest to point p.
		 */
		GEO2::index_t nearest_facet( const GEO2::vec3& p, GEO2::vec3& nearest_point, double& sq_dist ) const
		{
			GEO2::index_t nearest_facet;
			get_nearest_facet_hint( p, nearest_facet, nearest_point, sq_dist );
			nearest_facet_recursive( p, nearest_facet, nearest_point, sq_dist, 1, 0, mesh_.countTriangles() );
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
			{
				get_nearest_facet_hint( p, nearest_facet, nearest_point, sq_dist );
			}
			nearest_facet_recursive( p, nearest_facet, nearest_point, sq_dist, 1, 0, mesh_.countTriangles() );
		}

		/*
		 * Finds the nearest facet on the surface, but stops early if a
		 * point within a given distance is found.
		 */
		GEO2::index_t facet_in_envelope( const GEO2::vec3& p, double sq_epsilon, GEO2::vec3& nearest_point, double& sq_dist ) const
		{
			GEO2::index_t nearest_facet;
			get_nearest_facet_hint( p, nearest_facet, nearest_point, sq_dist );
			__m128i vec = _mm_setr_epi32( 1, 0, (int)mesh_.countTriangles(), 0 );
			__m256d pt = AvxMath::loadDouble3( &p.x );
			facet_in_envelope_recursive( pt, sq_epsilon, nearest_facet, nearest_point, sq_dist, vec );
			return nearest_facet;
		}

		/*
		 * Same as before, but stops as soon as a point on the surface in
		 * within a given distance bound from the triangle mesh.
		 */
		void facet_in_envelope_with_hint( const GEO2::vec3& p, double sq_epsilon, GEO2::index_t& nearest_facet, GEO2::vec3& nearest_point, double& sq_dist ) const
		{
			if( nearest_facet == GEO2::NO_FACET )
				get_nearest_facet_hint( p, nearest_facet, nearest_point, sq_dist );
			__m128i vec = _mm_setr_epi32( 1, 0, (int)mesh_.countTriangles(), 0 );
			__m256d pt = AvxMath::loadDouble3( &p.x );
			facet_in_envelope_recursive( pt, sq_epsilon, nearest_facet, nearest_point, sq_dist, vec );
		}

	  protected:
		/**
		 * \brief Computes a reasonable initialization for
		 *  nearest facet search.
		 *
		 * \details A good initialization makes the algorithm faster,
		 *  by allowing early pruning of subtrees that provably
		 *  do not contain the nearest neighbor.
		 *
		 * \param[in] p query point
		 * \param[out] nearest_facet a facet reasonably near p
		 * \param[out] nearest_point a point in nearest_facet
		 * \param[out] sq_dist squared distance between p and nearest_point
		 */
		void get_nearest_facet_hint( const GEO2::vec3& p, GEO2::index_t& nearest_facet, GEO2::vec3& nearest_point, double& sq_dist ) const;

		/**
		 * \brief The recursive function used by the implementation
		 *  of nearest_facet().
		 *
		 * \details The first call may use get_nearest_facet_hint()
		 * to initialize nearest_facet, nearest_point and sq_dist,
		 * as done in nearest_facet().
		 *
		 * \param[in] p query point
		 * \param[in,out] nearest_facet the nearest facet so far,
		 * \param[in,out] nearest_point a point in nearest_facet
		 * \param[in,out] sq_dist squared distance between p and nearest_point
		 * \param[in] n index of the current node in the AABB tree
		 * \param[in] b index of the first facet in the subtree under node \p n
		 * \param[in] e one position past the index of the last facet in the
		 *  subtree under node \p n
		 */
		void nearest_facet_recursive(
		  const GEO2::vec3& p, GEO2::index_t& nearest_facet, GEO2::vec3& nearest_point, double& sq_dist, GEO2::index_t n, GEO2::index_t b, GEO2::index_t e ) const;

		/*
		 * Same as before, but stops early if a point within a given distance
		 * is found.
		 */
		void facet_in_envelope_recursive(
		  __m256d p, double sq_epsilon, GEO2::index_t& nearest_facet, GEO2::vec3& nearest_point, double& sq_dist, __m128i nbe ) const;

	  protected:
		std::vector<GEO2::Box> bboxes_;
		const GEO2::Mesh& mesh_;
	};

	inline void get_point_facet_nearest_point( const GEO2::Mesh& M, __m256d p, GEO2::index_t f, GEO2::vec3& nearest_p, double& squared_dist )
	{
		using namespace GEO2;
		const vec3 *p1, *p2, *p3;
		M.getTriangleVertices( f, &p1, &p2, &p3 );
		squared_dist = point_triangle_squared_distance( p, *p1, *p2, *p3, &nearest_p );
	}
}  // namespace floatTetWild