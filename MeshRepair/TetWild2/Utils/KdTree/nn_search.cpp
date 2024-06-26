// Modified by Konstantin in late 2022, const.me, to make it stand-alone
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

#include <stdafx.h>
#include "nn_search.h"

namespace GEO2
{
	NearestNeighborSearch::NearestNeighborSearch( coord_index_t dimension )
		: dimension_( dimension )
		, nb_points_( 0 )
		, stride_( 0 )
		, points_( nullptr )
		, exact_( true )
	{
	}

	void NearestNeighborSearch::get_nearest_neighbors(
	  index_t nb_neighbors, const double* query_point, index_t* neighbors, double* neighbors_sq_dist, KeepInitialValues ) const
	{
		get_nearest_neighbors( nb_neighbors, query_point, neighbors, neighbors_sq_dist );
	}

	void NearestNeighborSearch::get_nearest_neighbors( index_t nb_neighbors, index_t query_point, index_t* neighbors, double* neighbors_sq_dist ) const
	{
		get_nearest_neighbors( nb_neighbors, point_ptr( query_point ), neighbors, neighbors_sq_dist );
	}

	void NearestNeighborSearch::set_points( index_t nb_points, const double* points )
	{
		nb_points_ = nb_points;
		points_ = points;
		stride_ = dimension_;
	}

	bool NearestNeighborSearch::stride_supported() const
	{
		return false;
	}

	void NearestNeighborSearch::set_points( index_t nb_points, const double* points, index_t stride )
	{
		if( stride == index_t( dimension() ) )
		{
			set_points( nb_points, points );
			return;
		}
		assert( stride_supported() );
		nb_points_ = nb_points;
		points_ = points;
		stride_ = stride;
	}

	void NearestNeighborSearch::set_exact( bool x )
	{
		exact_ = x;
	}
}  // namespace GEO2