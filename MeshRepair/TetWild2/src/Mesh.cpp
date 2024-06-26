// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "stdafx.h"
#include "Mesh.h"
#include "LocalOperations.h"
#include <numeric>
#include <omp.h>

namespace floatTetWild
{
	InsertionBuffers& Mesh::insertionBuffers() const
	{
#if PARALLEL_TRIANGLES_INSERTION
		return insertion[ omp_get_thread_num() ];
#else
		return insertion;
#endif
	}

	CollapseEdgeBuffers& Mesh::collapseEdgeBuffers()
	{
#if PARALLEL_TRIANGLES_INSERTION
		return collapseEdge[ omp_get_thread_num() ];
#else
		return collapseEdge;
#endif
	}

	void Mesh::one_ring_vertex_coloring( std::vector<Scalar>& colors ) const
	{
		colors.resize( tet_vertices.size() );
		std::fill( colors.begin(), colors.end(), -1 );
		colors[ 0 ] = 0;

		std::vector<bool> available( tet_vertices.size(), true );

		std::vector<int> ring;

		for( int i = 1; i < tet_vertices.size(); ++i )
		{
			const auto& v = tet_vertices[ i ];
			if( v.isRemoved() )
				continue;

			ring.clear();
			for( const auto& t : v.connTets )
			{
				for( int j = 0; j < 4; ++j )
					ring.push_back( tets[ t ][ j ] );
			}
			vector_unique( ring );

			for( const auto n : ring )
			{
				if( colors[ n ] != -1 )
					available[ colors[ n ] ] = false;
			}

			int first_available_col;
			for( first_available_col = 0; first_available_col < available.size(); first_available_col++ )
			{
				if( available[ first_available_col ] )
					break;
			}

			assert( available[ first_available_col ] );

			colors[ i ] = first_available_col;

			for( const auto n : ring )
			{
				if( colors[ n ] != -1 )
					available[ colors[ n ] ] = true;
			}
		}
	}

	void Mesh::one_ring_vertex_sets( const int threshold, std::vector<std::vector<int>>& concurrent_sets, std::vector<int>& serial_set ) const
	{
		std::vector<Scalar> colors;
		one_ring_vertex_coloring( colors );
		int max_c = -1;
		for( const auto c : colors )
			max_c = std::max( max_c, int( c ) );

		concurrent_sets.clear();
		concurrent_sets.resize( max_c + 1 );
		serial_set.clear();

		for( size_t i = 0; i < colors.size(); ++i )
		{
			const int col = colors[ i ];
			// removed vertex
			if( col < 0 )
				serial_set.push_back( i );
			else
				concurrent_sets[ col ].push_back( i );
		}

		for( int i = concurrent_sets.size() - 1; i >= 0; --i )
		{
			if( concurrent_sets[ i ].size() < threshold )
			{
				serial_set.insert( serial_set.end(), concurrent_sets[ i ].begin(), concurrent_sets[ i ].end() );
				concurrent_sets.erase( concurrent_sets.begin() + i );
			}
		}
	}

	void Mesh::one_ring_edge_set( const std::vector<std::array<int, 2>>& edges, const std::vector<bool>& v_is_removed, const std::vector<bool>& f_is_removed,
	  const std::vector<std::unordered_set<int>>& conn_fs, const std::vector<Vector3>& input_vertices, std::vector<int>& safe_set )
	{
		// std::vector<int> indices(edges.size());
		// std::vector<Scalar> weights(edges.size());

		// for(size_t e_id = 0; e_id < edges.size(); ++e_id)
		// {
		// 	indices[e_id] = e_id;
		// 	weights[e_id] = (input_vertices[edges[e_id][0]] - input_vertices[edges[e_id][1]]).squaredNorm();
		// }

		// std::sort(indices.begin(), indices.end(), [&](const int i, const int j){ return weights[i] < weights[j]; });
		std::vector<int> indices( edges.size() );
		std::iota( std::begin( indices ), std::end( indices ), 0 );
		// floatTetWild::Random::shuffle( indices );

		std::vector<bool> unsafe_face( f_is_removed.size(), false );
		safe_set.clear();
		for( const int e_id : indices )
		{
			const auto e = edges[ e_id ];
			if( v_is_removed[ e[ 0 ] ] || v_is_removed[ e[ 1 ] ] )
				continue;

			bool ok = true;

			for( const int f : conn_fs[ e[ 0 ] ] )
			{
				if( f_is_removed[ f ] )
					continue;

				if( unsafe_face[ f ] )
				{
					ok = false;
					break;
				}
			}
			if( !ok )
				continue;
			for( const int f : conn_fs[ e[ 1 ] ] )
			{
				if( f_is_removed[ f ] )
					continue;

				if( unsafe_face[ f ] )
				{
					ok = false;
					break;
				}
			}
			if( !ok )
				continue;

			safe_set.push_back( e_id );

			for( const int f : conn_fs[ e[ 0 ] ] )
			{
				if( f_is_removed[ f ] )
					continue;

				assert( !unsafe_face[ f ] );
				unsafe_face[ f ] = true;
			}
			for( const int f : conn_fs[ e[ 1 ] ] )
			{
				if( f_is_removed[ f ] )
					continue;

				// assert(!unsafe_face[f]);
				unsafe_face[ f ] = true;
			}
		}
	}

	void Mesh::createThreadLocalBuffers()
	{
		const size_t count = std::max( params.num_threads, 1u );
		findNewPosBuffers.resize( count );
		facetRecursionStacks.resize( count );
#if PARALLEL_TRIANGLES_INSERTION
		insertion.resize( count );
		collapseEdge.resize( count );
#endif
	}
}  // namespace floatTetWild