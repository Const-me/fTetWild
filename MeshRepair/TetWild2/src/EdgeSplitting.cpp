// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
#include "stdafx.h"
#include "EdgeSplitting.h"
#include "LocalOperations.h"
#include "MeshImprovement.h"

constexpr uint8_t TET_MODIFIED = 100;

void floatTetWild::edge_splitting( Mesh& mesh, const AABBWrapper& tree )
{
	auto tm = mesh.times.edgeSplitting.measure();
	auto& tets = mesh.tets;
	auto& tet_vertices = mesh.tet_vertices;

	int counter = 0;
	int suc_counter = 0;

	////init
	mesh.reset_t_empty_start();
	mesh.reset_v_empty_start();

	EdgesSet edges;
	get_all_edges( mesh, edges );

	std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_l> es_queue;

	edges.enumerate(
	  [ & ]( int e0, int e1 )
	  {
		  Scalar l_2 = get_edge_length_2( mesh, e0, e1 );
		  Scalar sizing_scalar = ( tet_vertices[ e0 ].sizing_scalar + tet_vertices[ e1 ].sizing_scalar ) / 2;
		  if( l_2 > mesh.params.split_threshold_2 * sizing_scalar * sizing_scalar )
			  es_queue.push( ElementInQueue( e0, e1, l_2 ) );
	  } );
	edges.clear();

	////split
	int budget = -1;  // input
	if( budget > 0 )
	{
		int v_slots = mesh.v_empty_size();
		v_slots = budget - v_slots;
		if( v_slots > 0 )
		{
			tet_vertices.reserve( tet_vertices.size() + v_slots );
		}
		int t_slots = mesh.t_empty_size();
		t_slots = budget * 6 - t_slots;
		if( t_slots > 0 )
		{
			tets.reserve( tet_vertices.size() + t_slots );
		}
	}
	else
	{
		// reserve space
		int v_slots = mesh.v_empty_size();
		int t_slots = mesh.t_empty_size();
		if( v_slots < es_queue.size() * 2 )
			tet_vertices.reserve( tet_vertices.size() + es_queue.size() * 2 - v_slots );
		if( t_slots < es_queue.size() * 6 * 2 )
			tets.reserve( tets.size() + es_queue.size() * 6 * 2 - t_slots + 1 );
	}

	std::vector<bool> is_splittable( mesh.tets.size(), true );
	bool is_repush = true;
	while( !es_queue.empty() )
	{
		std::array<int, 2> v_ids = es_queue.top().v_ids;
		es_queue.pop();

		if( tet_vertices[ v_ids[ 0 ] ].isFreezed() && tet_vertices[ v_ids[ 1 ] ].isFreezed() )
			continue;

		std::vector<std::array<int, 2>> new_edges;
		if( split_an_edge( mesh, v_ids[ 0 ], v_ids[ 1 ], is_repush, new_edges, is_splittable, tree ) )
			suc_counter++;
		else
			is_repush = false;

		for( auto& e : new_edges )
		{
			Scalar l_2 = get_edge_length_2( mesh, e[ 0 ], e[ 1 ] );
			Scalar sizing_scalar = ( tet_vertices[ e[ 0 ] ].sizing_scalar + tet_vertices[ e[ 1 ] ].sizing_scalar ) / 2;
			if( l_2 > mesh.params.split_threshold_2 * sizing_scalar * sizing_scalar )
			{
				es_queue.push( ElementInQueue( e, l_2 ) );
			}
		}

		counter++;

		if( budget > 0 )
		{
			budget--;
			if( budget == 0 )
				break;
		}
	}

	for( auto& t : tets )
	{
		if( t.is_removed )
			continue;
		if( t.scalar == TET_MODIFIED )
		{
			t.quality = get_quality( mesh, t );
			t.scalar = 0;
		}
	}

	mesh.logger().logDebug( "success = %i ( %i )", suc_counter, counter );
}

bool floatTetWild::split_an_edge(
  Mesh& mesh, int v1_id, int v2_id, bool is_repush, std::vector<std::array<int, 2>>& new_edges, std::vector<bool>& is_splittable, const AABBWrapper& tree )
{
	auto& tet_vertices = mesh.tet_vertices;
	auto& tets = mesh.tets;

	////create new vertex
	MeshVertex new_v;
	new_v.pos = ( tet_vertices[ v1_id ].pos + tet_vertices[ v2_id ].pos ) / 2;
	bool is_found = false;
	for( int i = mesh.v_empty_start; i < tet_vertices.size(); i++ )
	{
		mesh.v_empty_start = i;
		if( tet_vertices[ i ].isRemoved() )
		{
			is_found = true;
			break;
		}
	}
	if( !is_found )
		mesh.v_empty_start = tet_vertices.size();

	int v_id = mesh.v_empty_start;
	if( v_id < tet_vertices.size() )
		tet_vertices[ v_id ] = new_v;
	else
		tet_vertices.push_back( new_v );

	////check inversion
	std::vector<int> old_t_ids;
	setIntersection( tet_vertices[ v1_id ].connTets, tet_vertices[ v2_id ].connTets, old_t_ids );
	for( int t_id : old_t_ids )
	{
		if( !is_splittable[ t_id ] )
		{
			tet_vertices[ v_id ].setFlag( eVertexFlags::Removed );
			return false;
		}
		for( int j = 0; j < 4; j++ )
		{
			if( tets[ t_id ][ j ] == v1_id || tets[ t_id ][ j ] == v2_id )
			{
				if( is_inverted( mesh, t_id, j, new_v.pos ) )
				{
					for( int t_id1 : old_t_ids )
						is_splittable[ t_id1 ] = false;
					tet_vertices[ v_id ].setFlag( eVertexFlags::Removed );
					return false;
				}
			}
		}
	}

	////real update
	// update vertex
	tet_vertices[ v_id ].sizing_scalar = ( tet_vertices[ v1_id ].sizing_scalar + tet_vertices[ v2_id ].sizing_scalar ) / 2;
	tet_vertices[ v_id ].setFlag( eVertexFlags::BoundingBox, is_bbox_edge( mesh, v1_id, v2_id, old_t_ids ) );
	tet_vertices[ v_id ].setFlag( eVertexFlags::Surface, is_surface_edge( mesh, v1_id, v2_id, old_t_ids ) );
	tet_vertices[ v_id ].setFlag( eVertexFlags::Boundary, is_boundary_edge( mesh, v1_id, v2_id, tree ) );
	if( !mesh.is_input_all_inserted && tet_vertices[ v_id ].isBoundary() )
		tet_vertices[ v_id ].setFlag( eVertexFlags::Cut, tet_vertices[ v1_id ].isCut() && tet_vertices[ v2_id ].isCut() );

	// update tets
	std::vector<int> new_t_ids;
	get_new_tet_slots( mesh, old_t_ids.size(), new_t_ids );
	for( int t_id : new_t_ids )
		tets[ t_id ].reset();

	// update indices & tags
	for( int i = 0; i < old_t_ids.size(); i++ )
	{
		tets[ new_t_ids[ i ] ] = tets[ old_t_ids[ i ] ];
		for( int j = 0; j < 4; j++ )
		{
			if( tets[ old_t_ids[ i ] ][ j ] == v1_id )
				tets[ old_t_ids[ i ] ][ j ] = v_id;
			else if( tets[ old_t_ids[ i ] ][ j ] == v2_id )
			{
				tets[ old_t_ids[ i ] ].is_surface_fs[ j ] = NOT_SURFACE;
				tets[ old_t_ids[ i ] ].surface_tags[ j ] = NO_SURFACE_TAG;
				tets[ old_t_ids[ i ] ].is_bbox_fs[ j ] = NOT_BBOX;
			}

			if( tets[ new_t_ids[ i ] ][ j ] == v2_id )
				tets[ new_t_ids[ i ] ][ j ] = v_id;
			else if( tets[ new_t_ids[ i ] ][ j ] == v1_id )
			{
				tets[ new_t_ids[ i ] ].is_surface_fs[ j ] = NOT_SURFACE;
				tets[ new_t_ids[ i ] ].surface_tags[ j ] = NO_SURFACE_TAG;
				tets[ new_t_ids[ i ] ].is_bbox_fs[ j ] = NOT_BBOX;
			}
		}
	}
	// update quality
	for( int i = 0; i < old_t_ids.size(); i++ )
	{
		tets[ old_t_ids[ i ] ].scalar = TET_MODIFIED;
		tets[ new_t_ids[ i ] ].scalar = TET_MODIFIED;
	}

	tet_vertices[ v_id ].connTets.addRange( old_t_ids );
	tet_vertices[ v_id ].connTets.addRange( new_t_ids );
	for( int i = 0; i < old_t_ids.size(); i++ )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( tets[ old_t_ids[ i ] ][ j ] != v_id && tets[ old_t_ids[ i ] ][ j ] != v2_id )
				tet_vertices[ tets[ old_t_ids[ i ] ][ j ] ].connTets.add( new_t_ids[ i ] );
		}
		tet_vertices[ v1_id ].connTets.remove( old_t_ids[ i ] );
		tet_vertices[ v1_id ].connTets.add( new_t_ids[ i ] );
	}

	if( mesh.tets.size() != is_splittable.size() )
		is_splittable.resize( mesh.tets.size(), true );

	////repush
	if( !is_repush )
		return true;

	std::vector<int> n_v_ids = { v1_id, v2_id };
	n_v_ids.reserve( old_t_ids.size() * 3 );
	for( int t_id : old_t_ids )
	{
		for( int j = 0; j < 4; j++ )
			if( tets[ t_id ][ j ] != v_id )
				n_v_ids.push_back( tets[ t_id ][ j ] );
	}
	vector_unique( n_v_ids );

	for( int n_v_id : n_v_ids )
		new_edges.push_back( { { v_id, n_v_id } } );

	return true;
}