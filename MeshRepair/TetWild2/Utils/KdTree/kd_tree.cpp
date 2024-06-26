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
#include "kd_tree.h"
#include "../Geogram2.h"
#include "../../parallelThreadsCount.h"

namespace
{
	using namespace GEO2;
	using index_t = uint32_t;
	using coord_index_t = uint8_t;

	inline double geo_sqr( double x )
	{
		return x * x;
	}

	constexpr uint32_t max_index_t()
	{
		return UINT_MAX;
	}

	inline double distance2( const double* a, const double* b, uint8_t dims )
	{
		if( 3 == dims )
		{
			const vec3* av = (const vec3*)a;
			const vec3* bv = (const vec3*)b;
			return distance2( *av, *bv );
		}
		__debugbreak();
		return std::numeric_limits<double>::max();
	}

	/**
	 * \brief Comparison functor used to
	 *  sort the point indices. Used by
	 *  BalancedKdTree.
	 */
	class ComparePointCoord
	{
	  public:
		/**
		 * \brief Creates a new ComparePointCoord
		 * \param[in] nb_points number of points
		 * \param[in] points pointer to first point
		 * \param[in] stride number of doubles between two
		 *  consecutive points in array (=dimension if point
		 *  array is compact).
		 * \param[in] splitting_coord the coordinate to compare
		 */
		ComparePointCoord( index_t nb_points, const double* points, index_t stride, coord_index_t splitting_coord )
			: nb_points_( nb_points )
			, points_( points )
			, stride_( stride )
			, splitting_coord_( splitting_coord )
		{
		}

		/**
		 * \brief Compares to point indices (does the
		 * indirection and coordinate lookup).
		 * \param[in] i index of first point to compare
		 * \param[in] j index of second point to compare
		 * \return true if point \p i is before point \p j, false otherwise
		 */
		bool operator()( index_t i, index_t j ) const
		{
			assert( i < nb_points_ );
			assert( j < nb_points_ );
			return ( points_ + i * stride_ )[ splitting_coord_ ] < ( points_ + j * stride_ )[ splitting_coord_ ];
		}

	  private:
		index_t nb_points_;
		const double* points_;
		index_t stride_;
		coord_index_t splitting_coord_;
	};
}  // namespace

/****************************************************************************/

namespace GEO2
{
	KdTree::KdTree( coord_index_t dim )
		: NearestNeighborSearch( dim )
		, bbox_min_( dim )
		, bbox_max_( dim )
		, root_( index_t( -1 ) )
	{
	}

	bool KdTree::stride_supported() const
	{
		return true;
	}

	void KdTree::set_points( index_t nb_points, const double* points, index_t stride )
	{
		nb_points_ = nb_points;
		points_ = points;
		stride_ = stride;

		point_index_.resize( nb_points );
		for( index_t i = 0; i < nb_points; i++ )
		{
			point_index_[ i ] = i;
		}

		// Compute the bounding box.
		for( coord_index_t c = 0; c < dimension(); ++c )
		{
			bbox_min_[ c ] = std::numeric_limits<double>::max();
			bbox_max_[ c ] = std::numeric_limits<double>::lowest();
		}
		for( index_t i = 0; i < nb_points; ++i )
		{
			const double* p = point_ptr( i );
			for( coord_index_t c = 0; c < dimension(); ++c )
			{
				bbox_min_[ c ] = std::min( bbox_min_[ c ], p[ c ] );
				bbox_max_[ c ] = std::max( bbox_max_[ c ], p[ c ] );
			}
		}

		root_ = build_tree();
	}

	void KdTree::set_points( index_t nb_points, const double* points )
	{
		set_points( nb_points, points, dimension() );
	}

	void KdTree::get_nearest_neighbors( index_t nb_neighbors, const double* query_point, index_t* neighbors, double* neighbors_sq_dist ) const
	{
		assert( nb_neighbors <= nb_points() );

		// Compute distance between query point and global bounding box
		// and copy global bounding box to local variables (bbox_min, bbox_max),
		// allocated on the stack. bbox_min and bbox_max are updated during the
		// traversal of the BalancedKdTree (see
		// get_nearest_neighbors_recursive()). They are necessary to
		// compute the distance between the query point and the
		// bbox of the current node.
		double box_dist = 0.0;
		double* bbox_min = (double*)( alloca( dimension() * sizeof( double ) ) );
		double* bbox_max = (double*)( alloca( dimension() * sizeof( double ) ) );
		init_bbox_and_bbox_dist_for_traversal( bbox_min, bbox_max, box_dist, query_point );
		NearestNeighbors NN( nb_neighbors, neighbors, neighbors_sq_dist, (index_t*)alloca( sizeof( index_t ) * ( nb_neighbors + 1 ) ),
		  (double*)alloca( sizeof( double ) * ( nb_neighbors + 1 ) ) );
		get_nearest_neighbors_recursive( root_, 0, nb_points(), bbox_min, bbox_max, box_dist, query_point, NN );
		NN.copy_to_user();
	}

	void KdTree::get_nearest_neighbors(
	  index_t nb_neighbors, const double* query_point, index_t* neighbors, double* neighbors_sq_dist, KeepInitialValues KV ) const
	{
		assert( nb_neighbors <= nb_points() );
		// Compute distance between query point and global bounding box
		// and copy global bounding box to local variables (bbox_min, bbox_max),
		// allocated on the stack. bbox_min and bbox_max are updated during the
		// traversal of the BalancedKdTree
		// (see get_nearest_neighbors_recursive()). They
		// are necessary to compute the distance between the query point and the
		// bbox of the current node.
		double box_dist = 0.0;
		double* bbox_min = (double*)( alloca( dimension() * sizeof( double ) ) );
		double* bbox_max = (double*)( alloca( dimension() * sizeof( double ) ) );
		init_bbox_and_bbox_dist_for_traversal( bbox_min, bbox_max, box_dist, query_point );
		NearestNeighbors NN( nb_neighbors, neighbors, neighbors_sq_dist, (index_t*)alloca( sizeof( index_t ) * ( nb_neighbors + 1 ) ),
		  (double*)alloca( sizeof( double ) * ( nb_neighbors + 1 ) ) );
		NN.copy_from_user();
		get_nearest_neighbors_recursive( root_, 0, nb_points(), bbox_min, bbox_max, box_dist, query_point, NN );
		NN.copy_to_user();
	}

	void KdTree::get_nearest_neighbors( index_t nb_neighbors, index_t q_index, index_t* neighbors, double* neighbors_sq_dist ) const
	{
		// TODO: optimized version that uses the fact that
		// we know that query_point is in the search data
		// structure already.
		// (I tried something already, see in the Attic,
		//  but it did not give any significant speedup).
		get_nearest_neighbors( nb_neighbors, point_ptr( q_index ), neighbors, neighbors_sq_dist );
	}

	void KdTree::get_nearest_neighbors_recursive(
	  index_t node_index, index_t b, index_t e, double* bbox_min, double* bbox_max, double box_dist, const double* query_point, NearestNeighbors& NN ) const
	{
		assert( e > b );

		// Simple case (node is a leaf)
		if( ( e - b ) <= MAX_LEAF_SIZE )
		{
			get_nearest_neighbors_leaf( node_index, b, e, query_point, NN );
			return;
		}

		// Get node attributes (virtual function call).

		index_t left_node_index;
		index_t right_node_index;
		coord_index_t coord;
		index_t m;
		double val;

		get_node( node_index, b, e, left_node_index, right_node_index, coord, m, val );

		double cut_diff = query_point[ coord ] - val;

		// If the query point is on the left side
		if( cut_diff < 0.0 )
		{
			// Traverse left subtree
			{
				double bbox_max_save = bbox_max[ coord ];
				bbox_max[ coord ] = val;
				get_nearest_neighbors_recursive( left_node_index, b, m, bbox_min, bbox_max, box_dist, query_point, NN );
				bbox_max[ coord ] = bbox_max_save;
			}

			// Update bbox distance (now measures the
			// distance to the bbox of the right subtree)
			double box_diff = bbox_min[ coord ] - query_point[ coord ];
			if( box_diff > 0.0 )
			{
				box_dist -= geo_sqr( box_diff );
			}
			box_dist += geo_sqr( cut_diff );

			// Traverse the right subtree, only if bbox
			// distance is nearer than furthest neighbor,
			// else there is no chance that the right
			// subtree contains points that will change
			// anything in the nearest neighbors NN.
			if( box_dist <= NN.furthest_neighbor_sq_dist() )
			{
				double bbox_min_save = bbox_min[ coord ];
				bbox_min[ coord ] = val;
				get_nearest_neighbors_recursive( right_node_index, m, e, bbox_min, bbox_max, box_dist, query_point, NN );
				bbox_min[ coord ] = bbox_min_save;
			}
		}
		else
		{
			// else the query point is on the right side
			// (then do the same with left and right subtree
			//  permutted).
			{
				double bbox_min_save = bbox_min[ coord ];
				bbox_min[ coord ] = val;
				get_nearest_neighbors_recursive( right_node_index, m, e, bbox_min, bbox_max, box_dist, query_point, NN );
				bbox_min[ coord ] = bbox_min_save;
			}

			// Update bbox distance (now measures the
			// distance to the bbox of the left subtree)
			double box_diff = query_point[ coord ] - bbox_max[ coord ];
			if( box_diff > 0.0 )
			{
				box_dist -= geo_sqr( box_diff );
			}
			box_dist += geo_sqr( cut_diff );

			if( box_dist <= NN.furthest_neighbor_sq_dist() )
			{
				double bbox_max_save = bbox_max[ coord ];
				bbox_max[ coord ] = val;
				get_nearest_neighbors_recursive( left_node_index, b, m, bbox_min, bbox_max, box_dist, query_point, NN );
				bbox_max[ coord ] = bbox_max_save;
			}
		}
	}

	void KdTree::get_nearest_neighbors_leaf( index_t node_index, index_t b, index_t e, const double* query_point, NearestNeighbors& NN ) const
	{
		NN.nb_visited += ( e - b );
		double R = NN.furthest_neighbor_sq_dist();
		index_t nb = e - b;
		const index_t* idx = &point_index_[ b ];

		// TODO: check generated ASM (I'd like to have AVX here).
		// We may need to dispatch according to dimension.

		index_t local_idx[ MAX_LEAF_SIZE ];
		double local_sq_dist[ MAX_LEAF_SIZE ];

		// Cache indices and computed distances in local
		// array. I guess AVX likes that (to be checked).
		// Not sure, because access to p is indirect, maybe
		// I should copy the points to local memory before
		// computing the distances (or having another copy
		// of the points array that I pre-reorder so that
		// leaf's points are in a contiguous chunk of memory),
		// to be tested...
		for( index_t ii = 0; ii < nb; ++ii )
		{
			index_t i = idx[ ii ];
			const double* p = point_ptr( i );
			double sq_dist = ::distance2( query_point, p, dimension() );
			local_idx[ ii ] = i;
			local_sq_dist[ ii ] = sq_dist;
		}

		// Now insert the points that are nearer to query
		// point than NN's bounding ball.
		for( index_t ii = 0; ii < nb; ++ii )
		{
			double sq_dist = local_sq_dist[ ii ];
			if( sq_dist <= R )
			{
				NN.insert( local_idx[ ii ], sq_dist );
				R = NN.furthest_neighbor_sq_dist();
			}
		}
	}

	void KdTree::init_bbox_and_bbox_dist_for_traversal( double* bbox_min, double* bbox_max, double& box_dist, const double* query_point ) const
	{
		// Compute distance between query point and global bounding box
		// and copy global bounding box to local variables (bbox_min, bbox_max),
		// allocated on the stack. bbox_min and bbox_max are updated during the
		// traversal of the KdTree (see get_nearest_neighbors_recursive()). They
		// are necessary to compute the distance between the query point and the
		// bbox of the current node.
		box_dist = 0.0;
		for( coord_index_t c = 0; c < dimension(); ++c )
		{
			bbox_min[ c ] = bbox_min_[ c ];
			bbox_max[ c ] = bbox_max_[ c ];
			if( query_point[ c ] < bbox_min_[ c ] )
			{
				box_dist += geo_sqr( bbox_min_[ c ] - query_point[ c ] );
			}
			else if( query_point[ c ] > bbox_max_[ c ] )
			{
				box_dist += geo_sqr( bbox_max_[ c ] - query_point[ c ] );
			}
		}
	}

	/****************************************************************************/

	BalancedKdTree::BalancedKdTree( coord_index_t dim )
		: KdTree( dim )
	{
	}

	void BalancedKdTree::build_tree_omp()
	{
		std::array<index_t, 9> arr;
		memset( arr.data(), 0xFF, sizeof( arr ) );

		// Create the first level of the tree
		arr[ 0 ] = 0;
		arr[ 8 ] = nb_points();
		arr[ 4 ] = split_kd_node( 1, arr[ 0 ], arr[ 8 ] );

		// Create the second level of the tree, using two threads
#pragma omp parallel for
		for( int i = 0; i < 2; i++ )
			arr[ i * 4 + 2 ] = split_kd_node( i + 2, arr[ i * 4 ], arr[ i * 4 + 4 ] );

		// Create the third level of the tree, using four threads
#pragma omp parallel for
		for( int i = 0; i < 4; i++ )
			arr[ i * 2 + 1 ] = split_kd_node( i + 4, arr[ i * 2 ], arr[ i * 2 + 2 ] );

		// Create the fourth level of the tree, using eight threads
#pragma omp parallel for
		for( int i = 0; i < 8; i++ )
			create_kd_tree_recursive( i + 8, arr[ i ], arr[ i + 1 ] );
	}

	index_t BalancedKdTree::build_tree()
	{
		index_t sz = max_node_index( 1, 0, nb_points() ) + 1;
		splitting_coord_.resize( sz );
		splitting_val_.resize( sz );

		// If there are more than 16*MAX_LEAF_SIZE (=256) points,
		// create the tree in parallel
		if( nb_points() >= ( 16 * MAX_LEAF_SIZE ) && MeshRepair::shouldUseOpenMP() )
			build_tree_omp();
		else
			create_kd_tree_recursive( 1, 0, nb_points() );

		// Root node is number 1.
		// This is because "children at 2*n and 2*n+1" does not
		// work with 0 !!
		return 1;
	}

	index_t BalancedKdTree::split_kd_node( index_t node_index, index_t b, index_t e )
	{
		assert( e > b );
		// Do not split leafs
		if( b + 1 == e )
		{
			return b;
		}

		coord_index_t splitting_coord = best_splitting_coord( b, e );
		index_t m = b + ( e - b ) / 2;
		assert( m < e );

		// sorts the indices in such a way that points's
		// coordinates splitting_coord in [b,m) are smaller
		// than m's and points in [m,e) are
		// greater or equal to m's
		std::nth_element( point_index_.begin() + std::ptrdiff_t( b ), point_index_.begin() + std::ptrdiff_t( m ), point_index_.begin() + std::ptrdiff_t( e ),
		  ComparePointCoord( nb_points_, points_, stride_, splitting_coord ) );

		// Initialize node's variables (splitting coord and
		// splitting value)
		splitting_coord_[ node_index ] = splitting_coord;
		splitting_val_[ node_index ] = point_ptr( point_index_[ m ] )[ splitting_coord ];
		return m;
	}

	coord_index_t BalancedKdTree::best_splitting_coord( index_t b, index_t e )
	{
		// Returns the coordinates that maximizes
		// point's spread. We should probably
		// use a tradeoff between spread and
		// bbox shape ratio, as done in ANN, but
		// this simple method seems to give good
		// results in our case.
		coord_index_t result = 0;
		double max_spread = spread( b, e, 0 );
		for( coord_index_t c = 1; c < dimension(); ++c )
		{
			double coord_spread = spread( b, e, c );
			if( coord_spread > max_spread )
			{
				result = c;
				max_spread = coord_spread;
			}
		}
		return result;
	}

	void BalancedKdTree::get_node(
	  index_t n, index_t b, index_t e, index_t& left_child, index_t& right_child, coord_index_t& splitting_coord, index_t& m, double& splitting_val ) const
	{
		left_child = 2 * n;
		right_child = 2 * n + 1;
		splitting_coord = splitting_coord_[ n ];
		m = b + ( e - b ) / 2;
		splitting_val = splitting_val_[ n ];
	}

	/**************************************************************************/

	AdaptiveKdTree::AdaptiveKdTree( coord_index_t dim )
		: KdTree( dim )
	{
	}

	index_t AdaptiveKdTree::new_node()
	{
		splitting_coord_.push_back( 0 );
		splitting_val_.push_back( 0.0 );
		node_m_.push_back( 0 );
		node_right_child_.push_back( 0 );
		return nb_nodes() - 1;
	}

	index_t AdaptiveKdTree::build_tree()
	{
		// Create kd-tree. Local copy of the bbox is used, because it
		// is modified during traversal.
		double* bbox_min = (double*)( alloca( dimension() * sizeof( double ) ) );
		double* bbox_max = (double*)( alloca( dimension() * sizeof( double ) ) );
		for( coord_index_t c = 0; c < dimension(); ++c )
		{
			bbox_min[ c ] = bbox_min_[ c ];
			bbox_max[ c ] = bbox_max_[ c ];
		}

		splitting_coord_.resize( 0 );
		splitting_val_.resize( 0 );
		node_m_.resize( 0 );
		node_right_child_.resize( 0 );

		return create_kd_tree_recursive( 0, nb_points(), bbox_min, bbox_max );
	}

	index_t AdaptiveKdTree::create_kd_tree_recursive( index_t b, index_t e, double* bbox_min, double* bbox_max )
	{
		if( e - b <= MAX_LEAF_SIZE )
		{
			return index_t( -1 );
		}

		index_t m;
		coord_index_t cut_dim;
		double cut_val;

		//   Compute m, cut_dim and cut_val,
		// and reorganize indices along cut_val.
		split_kd_node( b, e, bbox_min, bbox_max, m, cut_dim, cut_val );

		index_t n = new_node();
		splitting_coord_[ n ] = cut_dim;
		splitting_val_[ n ] = cut_val;
		node_m_[ n ] = m;

		{
			double bbox_max_save = bbox_max[ cut_dim ];
			bbox_max[ cut_dim ] = cut_val;
			// This creates the left child. It does not need to be stored
			// because the tree is created in an order such that the
			// left child of node n is n+1.
			create_kd_tree_recursive( b, m, bbox_min, bbox_max );
			bbox_max[ cut_dim ] = bbox_max_save;
		}

		{
			double bbox_min_save = bbox_min[ cut_dim ];
			bbox_min[ cut_dim ] = cut_val;
			// Note: right_child needs to be copied to local variables
			// before being set in node_right_child_[n] because
			// create_kd_tree_recursive() modifies node_right_child_
			// (reallocates). If the following two lines are
			// done in a single assignment, then the computed reference to
			// node_right_child_[n]
			// in the lhs is no longer valid after the rhs is evaluated !
			index_t right_child = create_kd_tree_recursive( m, e, bbox_min, bbox_max );
			node_right_child_[ n ] = right_child;
			bbox_min[ cut_dim ] = bbox_min_save;
		}

		return n;
	}

	void AdaptiveKdTree::split_kd_node( index_t b, index_t e, double* bbox_min, double* bbox_max, index_t& m, coord_index_t& cut_dim, double& cut_val )
	{
		// Like "sliding midpoint split" in ANN.

		const double ERR = 0.001;

		// Find length of longest box size
		double max_length = -1.0;
		for( coord_index_t d = 0; d < dimension(); ++d )
		{
			double length = bbox_max[ d ] - bbox_min[ d ];
			max_length = std::max( max_length, length );
		}

		// Cutting coordinate
		cut_dim = 0;

		// Find long side with most spread
		double max_spread = -1.0;
		for( coord_index_t d = 0; d < dimension(); ++d )
		{
			double length = bbox_max[ d ] - bbox_min[ d ];
			// Is it among longest ?
			if( length >= ( 1.0 - ERR ) * max_length )
			{
				double spr = spread( b, e, d );
				if( spr > max_spread )
				{
					max_spread = spr;
					cut_dim = d;
				}
			}
		}

		double ideal_cut_val = 0.5 * ( bbox_min[ cut_dim ] + bbox_max[ cut_dim ] );

		double coord_min, coord_max;
		get_minmax( b, e, cut_dim, coord_min, coord_max );

		cut_val = ideal_cut_val;

		// Make it slide if need be.
		if( ideal_cut_val < coord_min )
		{
			cut_val = coord_min;
		}
		else if( ideal_cut_val > coord_max )
		{
			cut_val = coord_max;
		}

		index_t br1, br2;
		plane_split( b, e, cut_dim, cut_val, br1, br2 );

		index_t m0 = b + ( e - b ) / 2;
		m = m0;

		if( ideal_cut_val < coord_min )
		{
			m = b + 1;
		}
		else if( ideal_cut_val > coord_max )
		{
			m = e - 1;
		}
		else if( br1 > m0 )
		{
			m = br1;
		}
		else if( br2 < m0 )
		{
			m = br2;
		}
	}

	void AdaptiveKdTree::plane_split( index_t b_in, index_t e_in, coord_index_t coord, double val, index_t& br1_out, index_t& br2_out )
	{
		int b = int( b_in );
		int e = int( e_in );
		int l = b;
		int r = e - 1;
		while( true )
		{
			while( l < e && point_coord( l, coord ) < val )
			{
				++l;
			}
			while( r >= 0 && point_coord( r, coord ) >= val )
			{
				--r;
			}
			if( l > r )
			{
				break;
			}
			std::swap( point_index_[ l ], point_index_[ r ] );
			++l;
			--r;
		}
		int br1 = l;
		r = e - 1;
		while( true )
		{
			while( l < e && point_coord( l, coord ) <= val )
			{
				++l;
			}
			while( r >= br1 && point_coord( r, coord ) > val )
			{
				--r;
			}
			if( l > r )
			{
				break;
			}
			std::swap( point_index_[ l ], point_index_[ r ] );
			++l;
			--r;
		}
		int br2 = l;
		br1_out = index_t( br1 );
		br2_out = index_t( br2 );
	}

	void AdaptiveKdTree::get_node(
	  index_t n, index_t b, index_t e, index_t& left_child, index_t& right_child, coord_index_t& splitting_coord, index_t& m, double& splitting_val ) const
	{
		assert( n < nb_nodes() );
		left_child = n + 1;
		right_child = node_right_child_[ n ];
		splitting_coord = splitting_coord_[ n ];
		m = node_m_[ n ];
		splitting_val = splitting_val_[ n ];
	}

	/*************************************************************************/

}  // namespace GEO2
