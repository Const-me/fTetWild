// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "stdafx.h"
#include "EdgeCollapsing.h"
#include "LocalOperations.h"
#include "MeshImprovement.h"
#include "EdgesSet.h"
#define EC_POSTPROCESS true

namespace floatTetWild
{
	enum struct eCollapseStatus : int
	{
		FailInversion = -1,
		FailQuality = -2,
		FailEnvelope0 = -3,
		failEnvelope1 = -4,
		failEnvelope2 = -5,
		failEnvelope3 = -6,
		Success = 1,
		SuccessEnvelope = 2,
	};

	bool isSuccessStatus( eCollapseStatus cs )
	{
		return (int)cs > 0;
	}

	namespace
	{
		void edge_collapsing_aux( Mesh& mesh, const AABBWrapper& tree, EdgesSet& edges )
		{
			auto& tets = mesh.tets;
			auto& tet_vertices = mesh.tet_vertices;

			int counter = 0;
			int suc_counter = 0;
			int suc_counter_env = 0;
			EdgeCollapsingAuxBuffers& buffers = mesh.edgeCollapsingAuxBuffers;

			////init
			std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_s>& ec_queue = buffers.ec_queue;

			edges.enumerate(
			  [ & ]( int e0, int e1 )
			  {
				  Scalar l_2 = get_edge_length_2( mesh, e0, e1 );
				  if( is_collapsable_length( mesh, e0, e1, l_2 ) && is_collapsable_boundary( mesh, e0, e1, tree ) )
				  {
					  ec_queue.push( ElementInQueue( e0, e1, l_2 ) );
					  ec_queue.push( ElementInQueue( e1, e0, l_2 ) );
				  }
			  } );
			edges.clear();

			////collapse
			int ts = 0;

			std::vector<std::array<int, 2>>& inf_es = buffers.inf_es;
			inf_es.clear();

			std::vector<int>& inf_e_tss = buffers.inf_e_tss;
			inf_e_tss.clear();

			std::vector<int>& tet_tss = buffers.tet_tss;
			tet_tss.resize( tets.size() );
			if( !tet_tss.empty() )
				memset( tet_tss.data(), 0, tet_tss.size() * 4 );

			do
			{
				counter = 0;
				suc_counter = 0;
				suc_counter_env = 0;
				while( !ec_queue.empty() )
				{
					std::array<int, 2> v_ids = ec_queue.top().v_ids;
					Scalar old_weight = ec_queue.top().weight;
					ec_queue.pop();

					while( !ec_queue.empty() )
					{
						if( ec_queue.top().v_ids == v_ids )
							ec_queue.pop();
						else
							break;
					}

					if( is_edge_freezed( mesh, v_ids[ 0 ], v_ids[ 1 ] ) )
						continue;

					if( !is_valid_edge( mesh, v_ids[ 0 ], v_ids[ 1 ] ) )
						continue;

					if( !is_collapsable_boundary( mesh, v_ids[ 0 ], v_ids[ 1 ], tree ) )
						continue;

					Scalar weight = get_edge_length_2( mesh, v_ids[ 0 ], v_ids[ 1 ] );
					if( weight != old_weight || !is_collapsable_length( mesh, v_ids[ 0 ], v_ids[ 1 ], weight ) )
						continue;

					if( !is_collapsable_bbox( mesh, v_ids[ 0 ], v_ids[ 1 ] ) )
						continue;

					EdgesSet& new_edges = buffers.new_edges;
					new_edges.clear();

					const eCollapseStatus result = collapse_an_edge( mesh, v_ids[ 0 ], v_ids[ 1 ], tree, new_edges, ts, tet_tss );
					if( isSuccessStatus( result ) )
					{
						if( result == eCollapseStatus::SuccessEnvelope )
							suc_counter_env++;
						suc_counter++;

						new_edges.enumerate(
						  [ & ]( int e0, int e1 )
						  {
							  if( is_edge_freezed( mesh, e0, e1 ) )
								  return;
							  Scalar l_2 = get_edge_length_2( mesh, e0, e1 );
							  if( is_collapsable_length( mesh, e0, e1, l_2 ) )
							  {
								  ec_queue.push( ElementInQueue( e0, e1, l_2 ) );
								  ec_queue.push( ElementInQueue( e1, e0, l_2 ) );
							  }
						  } );
					}
#if EC_POSTPROCESS
					else
					{
						inf_es.push_back( v_ids );
						inf_e_tss.push_back( ts );
					}
#endif
					counter++;
				}

				mesh.logger().logDebug( "success(env) = %i", suc_counter_env );
				mesh.logger().logDebug( "success = %i ( %i )", suc_counter, counter );

#if EC_POSTPROCESS
				if( suc_counter == 0 )
#endif
					break;

#if EC_POSTPROCESS
				////postprocess
				std::vector<std::array<int, 2>>& tmp_inf_es = buffers.tmp_inf_es;
				tmp_inf_es.clear();
				const size_t inf_es_size = inf_es.size();
				tmp_inf_es.reserve( inf_es_size / 4 + 1 );
				for( size_t i = 0; i < inf_es_size; i++ )
				{
					const std::array<int, 2> inf_es_i = inf_es[ i ];
					if( is_edge_freezed( mesh, inf_es_i[ 0 ], inf_es_i[ 1 ] ) )
						continue;
					if( !is_valid_edge( mesh, inf_es_i[ 0 ], inf_es_i[ 1 ] ) )
						continue;

					Scalar weight = get_edge_length_2( mesh, inf_es_i[ 0 ], inf_es_i[ 1 ] );
					if( !is_collapsable_length( mesh, inf_es_i[ 0 ], inf_es_i[ 1 ], weight ) )
						continue;

					if( !is_collapsable_bbox( mesh, inf_es_i[ 0 ], inf_es_i[ 1 ] ) )
						continue;

					bool is_recal = false;
					const auto& conn_tets = tet_vertices[ inf_es_i[ 0 ] ].connTets;
					const int recallLimit = inf_e_tss[ i ];
					for( int t_id : conn_tets )
					{
						if( tet_tss[ t_id ] > recallLimit )
						{
							is_recal = true;
							break;
						}
					}
					if( is_recal )
						ec_queue.push( ElementInQueue( inf_es_i, weight ) );
					else
						tmp_inf_es.push_back( inf_es_i );
				}
				vector_unique( tmp_inf_es );
				inf_es = tmp_inf_es;

				ts++;
				inf_e_tss.resize( inf_es.size() );
				std::fill( inf_e_tss.begin(), inf_e_tss.end(), ts );
#endif
			} while( suc_counter > 0 );
		}
	}  // namespace
}  // namespace floatTetWild

void floatTetWild::edge_collapsing( Mesh& mesh, const AABBWrapper& tree )
{
	auto tm = mesh.times.edgeCollapsing.measure();
	auto& tet_vertices = mesh.tet_vertices;
	auto& tets = mesh.tets;

	EdgesSet& edges = mesh.edgeCollapsingBuffers.edges;
	edges.clear();

	get_all_edges( mesh, edges );
	edge_collapsing_aux( mesh, tree, edges );
}

floatTetWild::eCollapseStatus floatTetWild::collapse_an_edge(
  Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree, EdgesSet& new_edges, int ts, std::vector<int>& tet_tss, bool is_check_quality, bool is_update_tss )
{
	auto tm = mesh.times.collapseAnEdge.measure();
	auto& tet_vertices = mesh.tet_vertices;
	auto& tets = mesh.tets;

	////check vertices
	// check isolate surface points
	if( tet_vertices[ v1_id ].isSurface() && is_isolate_surface_point( mesh, v1_id ) )
	{
		tet_vertices[ v1_id ].clearFlag( eVertexFlags::Surface | eVertexFlags::Boundary );
	}
	// check boundary/surface
	if( tet_vertices[ v1_id ].isBoundary() &&
		is_point_out_boundary_envelope( mesh, tet_vertices[ v2_id ].pos, tree ) )  // todo: you should check/unmark is_on_boundary around here
		return eCollapseStatus::FailEnvelope0;
	if( tet_vertices[ v1_id ].isSurface() && is_point_out_envelope( mesh, tet_vertices[ v2_id ].pos, tree ) )
		return eCollapseStatus::failEnvelope1;

	CollapseEdgeBuffers& buffers = mesh.collapseEdgeBuffers();
	////check tets
	std::vector<int>& n12_t_ids = buffers.n12_t_ids;
	n12_t_ids.clear();

	setIntersection( tet_vertices[ v1_id ].connTets, tet_vertices[ v2_id ].connTets, n12_t_ids );
	if( n12_t_ids.empty() )
		return eCollapseStatus::FailInversion;

	std::vector<int>& n1_t_ids = buffers.n1_t_ids;	// v1.conn_tets - n12_t_ids
	n1_t_ids.clear();

	tet_vertices[ v1_id ].connTets.sort();
	std::set_difference(
	  tet_vertices[ v1_id ].connTets.begin(), tet_vertices[ v1_id ].connTets.end(), n12_t_ids.begin(), n12_t_ids.end(), std::back_inserter( n1_t_ids ) );

	// inversion
	std::vector<int>& js_n1_t_ids = buffers.js_n1_t_ids;
	js_n1_t_ids.clear();

	for( int t_id : n1_t_ids )
	{
		int j = tets[ t_id ].find( v1_id );
		js_n1_t_ids.push_back( j );
		assert( j < 4 );
		if( is_inverted( mesh, t_id, j, tet_vertices[ v2_id ].pos ) )
			return eCollapseStatus::FailInversion;
	}

	// quality
	double old_max_quality = 0;
	if( mesh.is_coarsening )
		old_max_quality = mesh.params.stop_energy;
	else
	{
		if( is_check_quality )
		{
			for( int t_id : tet_vertices[ v1_id ].connTets )
			{
				if( tets[ t_id ].quality > old_max_quality )
					old_max_quality = tets[ t_id ].quality;
			}
		}
	}
	std::vector<Scalar>& new_qs = buffers.new_qs;
	new_qs.clear();
	new_qs.reserve( tet_vertices[ v1_id ].connTets.size() );
	int ii = 0;
	for( int t_id : n1_t_ids )
	{
		int j = js_n1_t_ids[ ii++ ];
		Scalar new_q = get_quality( tet_vertices[ v2_id ], tet_vertices[ tets[ t_id ][ mod4( j + 1 ) ] ], tet_vertices[ tets[ t_id ][ mod4( j + 2 ) ] ],
		  tet_vertices[ tets[ t_id ][ mod4( j + 3 ) ] ] );
		if( is_check_quality && new_q > old_max_quality )
			return eCollapseStatus::FailQuality;
		new_qs.push_back( new_q );
	}

	// envelope
	Scalar l = get_edge_length_2( mesh, v1_id, v2_id );
	if( l > 0 )
	{
		if( tet_vertices[ v1_id ].isBoundary() )
		{
			if( is_out_boundary_envelope( mesh, v1_id, tet_vertices[ v2_id ].pos, tree ) )
				return eCollapseStatus::failEnvelope2;
		}
		if( tet_vertices[ v1_id ].isSurface() )
		{
			if( is_out_envelope( mesh, v1_id, tet_vertices[ v2_id ].pos, tree ) )
				return eCollapseStatus::failEnvelope3;
		}
	}

	////real update
	// vertex
	tet_vertices[ v1_id ].setFlag( eVertexFlags::Removed );
	tet_vertices[ v2_id ].setFlag( eVertexFlags::BoundingBox, tet_vertices[ v1_id ].isBoundingBox() || tet_vertices[ v2_id ].isBoundingBox() );
	tet_vertices[ v2_id ].setFlag( eVertexFlags::Surface, tet_vertices[ v1_id ].isSurface() || tet_vertices[ v2_id ].isSurface() );
	tet_vertices[ v2_id ].setFlag( eVertexFlags::Boundary, tet_vertices[ v1_id ].isBoundary() || tet_vertices[ v2_id ].isBoundary() );

	// tets
	// update quality
	int i = 0;
	for( int t_id : n1_t_ids )
		tets[ t_id ].quality = new_qs[ i++ ];

	std::vector<int>& n1_v_ids = buffers.n1_v_ids;
	n1_v_ids.clear();
	n1_v_ids.reserve( n1_t_ids.size() * 4 );
	for( int t_id : n1_t_ids )
	{
		for( int j = 0; j < 4; j++ )
			n1_v_ids.push_back( tets[ t_id ][ j ] );
	}
	vector_unique( n1_v_ids );

	// update tags
	for( int t_id : n12_t_ids )
	{
		int8_t sf_facing_v1 = NOT_SURFACE;
		int8_t sf_facing_v2 = NOT_SURFACE;
		int tag_facing_v1 = NO_SURFACE_TAG;
		int tag_facing_v2 = NO_SURFACE_TAG;
		int bbox_facing_v1 = NOT_BBOX;
		int bbox_facing_v2 = NOT_BBOX;

		std::array<int, 2> j12;
		for( int j = 0; j < 4; j++ )
		{
			if( tets[ t_id ][ j ] == v1_id )
			{
				sf_facing_v1 = tets[ t_id ].is_surface_fs[ j ];
				tag_facing_v1 = tets[ t_id ].surface_tags[ j ];
				bbox_facing_v1 = tets[ t_id ].is_bbox_fs[ j ];
				j12[ 0 ] = j;
			}
			else if( tets[ t_id ][ j ] == v2_id )
			{
				sf_facing_v2 = tets[ t_id ].is_surface_fs[ j ];
				tag_facing_v2 = tets[ t_id ].surface_tags[ j ];
				bbox_facing_v2 = tets[ t_id ].is_bbox_fs[ j ];
				j12[ 1 ] = j;
			}
		}

		std::array<int, 2> sf_connecting_v12 = { { NOT_SURFACE, NOT_SURFACE } };
		std::array<int, 2> tag_connecting_v12 = { { NO_SURFACE_TAG, NO_SURFACE_TAG } };
		if( sf_facing_v2 != NOT_SURFACE && sf_facing_v1 != NOT_SURFACE )
		{
			//            sf_connecting_v12[0] = -sf_facing_v2 + sf_facing_v1;
			int new_tag = -sf_facing_v2 + sf_facing_v1;
			if( new_tag == 0 )
				//                sf_connecting_v12[0] = NOT_SURFACE;
				sf_connecting_v12[ 0 ] = 0;
			else if( new_tag > 0 )
				sf_connecting_v12[ 0 ] = 1;
			else
				sf_connecting_v12[ 0 ] = -1;
			tag_connecting_v12[ 0 ] = tag_facing_v1;
		}
		else if( sf_facing_v2 != NOT_SURFACE )
		{
			sf_connecting_v12[ 0 ] = -sf_facing_v2;
			tag_connecting_v12[ 0 ] = tag_facing_v2;
		}
		else if( sf_facing_v1 != NOT_SURFACE )
		{
			sf_connecting_v12[ 0 ] = sf_facing_v1;
			tag_connecting_v12[ 0 ] = tag_facing_v1;
		}
		if( sf_connecting_v12[ 0 ] != NOT_SURFACE )
		{
			sf_connecting_v12[ 1 ] = -sf_connecting_v12[ 0 ];
			tag_connecting_v12[ 1 ] = tag_connecting_v12[ 0 ];
		}

		std::array<int, 2> bbox_connecting_v12;
		bbox_connecting_v12[ 0 ] = NOT_BBOX;
		if( bbox_facing_v2 != NOT_BBOX )
			bbox_connecting_v12[ 0 ] = bbox_facing_v2;
		else if( sf_facing_v1 != NOT_BBOX )
			bbox_connecting_v12[ 0 ] = bbox_facing_v1;
		bbox_connecting_v12[ 1 ] = bbox_connecting_v12[ 0 ];

		for( int i = 0; i < 2; i++ )
		{
			std::vector<int>& pair = buffers.pair;
			pair.clear();
			setIntersection( tet_vertices[ tets[ t_id ][ mod4( j12[ mod2( i + 1 ) ] + 1 ) ] ].connTets,
			  tet_vertices[ tets[ t_id ][ mod4( j12[ mod2( i + 1 ) ] + 2 ) ] ].connTets,
			  tet_vertices[ tets[ t_id ][ mod4( j12[ mod2( i + 1 ) ] + 3 ) ] ].connTets, pair );
			if( pair.size() > 1 )
			{
				int opp_t_id = pair[ 0 ] == t_id ? pair[ 1 ] : pair[ 0 ];
				for( int j = 0; j < 4; j++ )
				{
					if( tets[ opp_t_id ][ j ] != tets[ t_id ][ mod4( j12[ mod2( i + 1 ) ] + 1 ) ] &&
						tets[ opp_t_id ][ j ] != tets[ t_id ][ mod4( j12[ mod2( i + 1 ) ] + 2 ) ] &&
						tets[ opp_t_id ][ j ] != tets[ t_id ][ mod4( j12[ mod2( i + 1 ) ] + 3 ) ] )
					{
						tets[ opp_t_id ].is_surface_fs[ j ] = sf_connecting_v12[ i ];
						tets[ opp_t_id ].surface_tags[ j ] = tag_connecting_v12[ i ];
						tets[ opp_t_id ].is_bbox_fs[ j ] = bbox_connecting_v12[ i ];
						break;
					}
				}
			}
		}
	}

	// update connectivity
	ts++;
	ii = 0;
	for( int t_id : n1_t_ids )
	{
		int j = js_n1_t_ids[ ii++ ];
		tets[ t_id ][ j ] = v2_id;
		tet_vertices[ v2_id ].connTets.add( t_id );
		if( is_update_tss )
			tet_tss[ t_id ] = ts;  // update timestamp
	}
	for( int t_id : n12_t_ids )
	{
		tets[ t_id ].is_removed = true;
		for( int j = 0; j < 4; j++ )
		{
			if( tets[ t_id ][ j ] != v1_id )
				tet_vertices[ tets[ t_id ][ j ] ].connTets.remove( t_id );
		}
	}

	tet_vertices[ v1_id ].connTets.clear();

	////re-push
	for( int v_id : n1_v_ids )
		if( v_id != v1_id )
			new_edges.add( v2_id, v_id );

	if( tet_vertices[ v1_id ].isSurface() )
		return eCollapseStatus::SuccessEnvelope;
	return eCollapseStatus::Success;
}

bool floatTetWild::is_edge_freezed( Mesh& mesh, int v1_id, int v2_id )
{
	if( mesh.tet_vertices[ v1_id ].isFreezed() || mesh.tet_vertices[ v2_id ].isFreezed() )
		return true;
	return false;
}

bool floatTetWild::is_collapsable_bbox( Mesh& mesh, int v1_id, int v2_id )
{
	if( !mesh.tet_vertices[ v1_id ].isBoundingBox() )
		return true;
	else if( !mesh.tet_vertices[ v2_id ].isBoundingBox() )
		return false;

	//    std::unordered_set<int> bbox_fs2;
	//    for (int t_id:mesh.tet_vertices[v2_id].conn_tets) {
	//        for (int j = 0; j < 4; j++) {
	//            if (mesh.tets[t_id][j] != v2_id && mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX)
	//                bbox_fs2.insert(mesh.tets[t_id].is_bbox_fs[j]);
	//        }
	//    }
	//    int old_size = bbox_fs2.size();
	//    for (int t_id:mesh.tet_vertices[v1_id].conn_tets) {
	//        for (int j = 0; j < 4; j++) {
	//            if (mesh.tets[t_id][j] != v1_id && mesh.tets[t_id].is_bbox_fs[j] != NOT_BBOX) {
	//                bbox_fs2.insert(mesh.tets[t_id].is_bbox_fs[j]);
	//                if (bbox_fs2.size() > old_size)
	//                    return false;
	//            }
	//        }
	//    }

	std::vector<int> bbox_fs2;
	for( int t_id : mesh.tet_vertices[ v2_id ].connTets )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( mesh.tets[ t_id ][ j ] != v2_id && mesh.tets[ t_id ].is_bbox_fs[ j ] != NOT_BBOX )
				bbox_fs2.push_back( mesh.tets[ t_id ].is_bbox_fs[ j ] );
		}
	}
	vector_unique( bbox_fs2 );

	for( int t_id : mesh.tet_vertices[ v1_id ].connTets )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( mesh.tets[ t_id ][ j ] != v1_id && mesh.tets[ t_id ].is_bbox_fs[ j ] != NOT_BBOX )
			{
				if( std::find( bbox_fs2.begin(), bbox_fs2.end(), mesh.tets[ t_id ].is_bbox_fs[ j ] ) == bbox_fs2.end() )
					return false;
			}
		}
	}

	return true;
}

bool floatTetWild::is_collapsable_length( Mesh& mesh, int v1_id, int v2_id, Scalar l_2 )
{
	Scalar sizing_scalar = ( mesh.tet_vertices[ v1_id ].sizing_scalar + mesh.tet_vertices[ v2_id ].sizing_scalar ) / 2;
	if( l_2 <= mesh.params.collapse_threshold_2 * sizing_scalar * sizing_scalar )
		return true;
	return false;
}

bool floatTetWild::is_collapsable_boundary( Mesh& mesh, int v1_id, int v2_id, const AABBWrapper& tree )
{
	if( mesh.tet_vertices[ v1_id ].isBoundary() && !is_boundary_edge( mesh, v1_id, v2_id, tree ) )
		return false;
	return true;

	//    if (mesh.tet_vertices[v1_id].on_boundary_e_id >= 0 && mesh.tet_vertices[v2_id].on_boundary_e_id
	//        && mesh.tet_vertices[v1_id].on_boundary_e_id != mesh.tet_vertices[v2_id].on_boundary_e_id)
	//        return false;
	//    return true;
}