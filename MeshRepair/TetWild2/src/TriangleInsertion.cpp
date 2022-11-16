// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

//
// Created by Yixin Hu on 2019-08-27.
//
#include "stdafx.h"
#include "TriangleInsertion.h"
#include "LocalOperations.h"
#include "../external/Predicates.h"
#include "CutTable2.h"
#include "intersections.h"
#include "MeshImprovement.h"  //fortest
#include <igl/Timer.h>
#include <Utils/randomShuffle.h>
#include <bitset>
#include <numeric>
#include <unordered_map>
#include "SmallBuffer.h"
#define III -1
#include "../external/Rational.h"
#include "../Utils/miscUtils.h"

floatTetWild::Vector3 floatTetWild::get_normal( const Vector3& a, const Vector3& b, const Vector3& c )
{
	return ( ( b - c ).cross( a - c ) ).normalized();
}

void floatTetWild::match_surface_fs( const Mesh& mesh, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  std::vector<bool>& is_face_inserted, std::vector<std::array<std::vector<int>, 4>>& track_surface_fs )
{
	auto comp = []( const std::array<int, 4>& a, const std::array<int, 4>& b )
	{ return std::tuple<int, int, int>( a[ 0 ], a[ 1 ], a[ 2 ] ) < std::tuple<int, int, int>( b[ 0 ], b[ 1 ], b[ 2 ] ); };

	std::vector<std::array<int, 4>> input_fs( input_faces.size() );
	for( int i = 0; i < input_faces.size(); i++ )
	{
		input_fs[ i ] = { { input_faces[ i ][ 0 ], input_faces[ i ][ 1 ], input_faces[ i ][ 2 ], i } };
		std::sort( input_fs[ i ].begin(), input_fs[ i ].begin() + 3 );
	}
	std::sort( input_fs.begin(), input_fs.end(), comp );

	for( int i = 0; i < mesh.tets.size(); i++ )
	{
		auto& t = mesh.tets[ i ];
		for( int j = 0; j < 4; j++ )
		{
			std::array<int, 3> f = { { t[ mod4( j + 1 ) ], t[ mod4( j + 2 ) ], t[ mod4( j + 3 ) ] } };
			std::sort( f.begin(), f.end() );
			auto bounds = std::equal_range( input_fs.begin(), input_fs.end(), std::array<int, 4>( { { f[ 0 ], f[ 1 ], f[ 2 ], -1 } } ), comp );
			for( auto it = bounds.first; it != bounds.second; ++it )
			{
				int f_id = ( *it )[ 3 ];
				is_face_inserted[ f_id ] = true;
				track_surface_fs[ i ][ j ].push_back( f_id );
			}
		}
	}
}

void floatTetWild::sort_input_faces(
  const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const Mesh& mesh, std::vector<int>& sorted_f_ids )
{
	std::vector<Scalar> weights( input_faces.size() );
	sorted_f_ids.resize( input_faces.size() );
	for( int i = 0; i < input_faces.size(); i++ )
	{
		sorted_f_ids[ i ] = i;
		Vector3 u = input_vertices[ input_faces[ i ][ 1 ] ] - input_vertices[ input_faces[ i ][ 0 ] ];
		Vector3 v = input_vertices[ input_faces[ i ][ 2 ] ] - input_vertices[ input_faces[ i ][ 0 ] ];
		weights[ i ] = u.cross( v ).squaredNorm();
	}

	if( mesh.params.not_sort_input )
		return;

	randomShuffle( sorted_f_ids );
}

void floatTetWild::insert_triangles( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::vector<int>& input_tags,
  Mesh& mesh, std::vector<bool>& is_face_inserted, AABBWrapper& tree, bool is_again )
{
	insert_triangles_aux( input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, is_again );
}

namespace
{
	using namespace floatTetWild;

	static void computeMeshQuality( Mesh& mesh )
	{
		for( auto& t : mesh.tets )
		{
			if( !t.is_removed )
				t.quality = get_quality( mesh, t );
		}
	}

	static void computeMeshQualityOmp( Mesh& mesh )
	{
		const int64_t length = (int64_t)mesh.tets.size();
#pragma omp parallel for schedule( dynamic )
		for( int64_t i = 0; i < length; i++ )
		{
			MeshTet& t = mesh.tets[ i ];
			if( !t.is_removed )
				t.quality = get_quality( mesh, t );
		}
	}
}  // namespace

void floatTetWild::insert_triangles_aux( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  const std::vector<int>& input_tags, Mesh& mesh, std::vector<bool>& is_face_inserted, AABBWrapper& tree, bool is_again )
{
	auto tm = mesh.times.insertTrianglesAux.measure();

	std::vector<bool> old_is_face_inserted = is_face_inserted;	/// is_face_inserted has been initialized in main

	mesh.logger().logInfo( "triangle insertion start, #f = %zu, #v = %zu, #t = %zu", input_faces.size(), mesh.tet_vertices.size(), mesh.tets.size() );
	/////
	std::vector<std::array<std::vector<int>, 4>> track_surface_fs( mesh.tets.size() );
	if( !is_again )
	{
		match_surface_fs( mesh, input_vertices, input_faces, is_face_inserted, track_surface_fs );
	}
	int cnt_matched = std::count( is_face_inserted.begin(), is_face_inserted.end(), true );
	mesh.logger().logInfo( "matched #f = %i, uninserted #f = %i", cnt_matched, (int)is_face_inserted.size() - cnt_matched );

	std::vector<int> sorted_f_ids;
	sort_input_faces( input_vertices, input_faces, mesh, sorted_f_ids );

	/////
	std::vector<Vector3> new_vertices;
	std::vector<std::array<int, 4>> new_tets;
	int cnt_fail = 0;
	int cnt_total = 0;

	//////
	TriangleInsertionVars& vars = mesh.globalVars.triangleInsertion;
	for( int i = 0; i < sorted_f_ids.size(); i++ )
	{
		// fortest
		if( !is_again && i > 0 && i % 1000 == 0 )
		{
			mesh.logger().logDebug( "inserting f%i... %i failed", i, cnt_fail );
			mesh.logger().logDebug( "snapped %i/%i", vars.cnt_snapped, cnt_total );
			mesh.logger().logDebug(
			  "\t- time_find_cutting_tets = %gs (total %gs)", vars.time_find_cutting_tets - vars.old_time_find_cutting_tets, vars.time_find_cutting_tets );
			mesh.logger().logDebug( "\t- time_cut_mesh = %gs (total %gs)", vars.time_cut_mesh - vars.old_time_cut_mesh, vars.time_cut_mesh );
			mesh.logger().logDebug( "\t- time_get_intersecting_edges_and_points = %gs (total %gs)",
			  vars.time_get_intersecting_edges_and_points - vars.old_time_get_intersecting_edges_and_points, vars.time_get_intersecting_edges_and_points );
			mesh.logger().logDebug(
			  "\t- time_subdivide_tets = %gs (total %gs)", vars.time_subdivide_tets - vars.old_time_subdivide_tets, vars.time_subdivide_tets );
			mesh.logger().logDebug(
			  "\t- time_push_new_tets = %gs (total %gs)", vars.time_push_new_tets - vars.old_time_push_new_tets, vars.time_push_new_tets );
			mesh.logger().logDebug( "\t- time_simplify_subdivision_result = %gs (total %gs)",
			  vars.time_simplify_subdivision_result - vars.old_time_simplify_subdivision_result, vars.time_simplify_subdivision_result );

			vars.old_time_find_cutting_tets = vars.time_find_cutting_tets;
			vars.old_time_cut_mesh = vars.time_cut_mesh;
			vars.old_time_get_intersecting_edges_and_points = vars.time_get_intersecting_edges_and_points;
			vars.old_time_subdivide_tets = vars.time_subdivide_tets;
			vars.old_time_push_new_tets = vars.time_push_new_tets;
			vars.old_time_simplify_subdivision_result = vars.time_simplify_subdivision_result;
			mesh.logger().logDebug( "#v = %i/%zu", mesh.get_v_num(), mesh.tet_vertices.size() );
			mesh.logger().logDebug( "#t = %i/%zu", mesh.get_t_num(), mesh.tets.size() );
		}

		int f_id = sorted_f_ids[ i ];
		if( is_face_inserted[ f_id ] )
			continue;

		cnt_total++;
		if( insert_one_triangle( f_id, input_vertices, input_faces, input_tags, mesh, track_surface_fs, tree, is_again ) )
			is_face_inserted[ f_id ] = true;
		else
			cnt_fail++;

		//        pausee();//fortest
		if( f_id == III )
			break;	// fortest
	}
	mesh.logger().logInfo( "insert_one_triangle * n done, #v = %zu, #t = %zu", mesh.tet_vertices.size(), mesh.tets.size() );
	mesh.logger().logInfo(
	  "uninserted #f = %zu/%i", std::count( is_face_inserted.begin(), is_face_inserted.end(), false ), (int)is_face_inserted.size() - cnt_matched );
	mesh.logger().logInfo( "total timing: %gs",
	  vars.time_find_cutting_tets + vars.time_cut_mesh + vars.time_get_intersecting_edges_and_points + vars.time_subdivide_tets + vars.time_push_new_tets );

	pair_track_surface_fs( mesh, track_surface_fs );
	mesh.logger().logInfo( "pair_track_surface_fs done" );

	/////
	std::vector<std::array<int, 2>> b_edges1;
	std::vector<std::pair<std::array<int, 2>, std::vector<int>>> b_edge_infos;
	std::vector<bool> is_on_cut_edges;
	find_boundary_edges( input_vertices, input_faces, is_face_inserted, old_is_face_inserted, b_edge_infos, is_on_cut_edges, b_edges1, mesh.logger() );
	mesh.logger().logInfo( "find_boundary_edges done" );
	std::vector<std::array<int, 3>> known_surface_fs;
	std::vector<std::array<int, 3>> known_not_surface_fs;
	insert_boundary_edges( input_vertices, input_faces, b_edge_infos, is_on_cut_edges, track_surface_fs, mesh, tree, is_face_inserted, is_again,
	  known_surface_fs, known_not_surface_fs );
	mesh.logger().logInfo(
	  "uninserted #f = %zu/%i", std::count( is_face_inserted.begin(), is_face_inserted.end(), false ), (int)is_face_inserted.size() - cnt_matched );

	std::vector<std::array<int, 2>> b_edges2;
	mark_surface_fs(
	  input_vertices, input_faces, input_tags, track_surface_fs, is_face_inserted, known_surface_fs, known_not_surface_fs, b_edges2, mesh, tree );
	mesh.logger().logInfo( "mark_surface_fs done" );

	// build b_tree using b_edges
	tree.init_tmp_b_mesh_and_tree( input_vertices, input_faces, b_edges1, mesh, b_edges2 );

	if( mesh.params.num_threads > 1 )
		computeMeshQualityOmp( mesh );
	else
		computeMeshQuality( mesh );

	if( std::count( is_face_inserted.begin(), is_face_inserted.end(), false ) == 0 )
		mesh.is_input_all_inserted = true;
	mesh.logger().logInfo( "#b_edge1 = %zu, #b_edges2 = %zu", b_edges1.size(), b_edges2.size() );
	pausee();
}

bool floatTetWild::insert_one_triangle( int insert_f_id, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  const std::vector<int>& input_tags, Mesh& mesh, std::vector<std::array<std::vector<int>, 4>>& track_surface_fs, AABBWrapper& tree, bool is_again )
{
	auto tm = mesh.times.insertOneTriangle.measure();

	std::array<Vector3, 3> vs = { { input_vertices[ input_faces[ insert_f_id ][ 0 ] ], input_vertices[ input_faces[ insert_f_id ][ 1 ] ],
	  input_vertices[ input_faces[ insert_f_id ][ 2 ] ] } };
	Vector3 n = ( vs[ 1 ] - vs[ 0 ] ).cross( vs[ 2 ] - vs[ 0 ] );
	n.normalize();
	int t = get_t( vs[ 0 ], vs[ 1 ], vs[ 2 ] );

	std::vector<int> cut_t_ids;
	find_cutting_tets( insert_f_id, input_vertices, input_faces, vs, mesh, cut_t_ids, is_again );

	if( cut_t_ids.empty() )
		throw std::logic_error( "cut_t_ids.empty()" );

	CutMesh cut_mesh( mesh, n, vs );
	cut_mesh.construct( cut_t_ids );

	if( cut_mesh.snap_to_plane() )
	{
		mesh.globalVars.triangleInsertion.cnt_snapped++;
		cut_mesh.project_to_plane( input_vertices.size() );
		cut_mesh.expand_new( cut_t_ids );
		cut_mesh.project_to_plane( input_vertices.size() );
	}

	std::vector<Vector3> points;
	std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
	std::vector<int> subdivide_t_ids;
	if( !cut_mesh.get_intersecting_edges_and_points( points, map_edge_to_intersecting_point, subdivide_t_ids ) )
	{
		if( is_again )
		{
			if( is_uninserted_face_covered( insert_f_id, input_vertices, input_faces, cut_t_ids, mesh ) )
				return true;
		}
		mesh.logger().logError( "FAIL get_intersecting_edges_and_points" );
		return false;
	}

	// have to add all cut_t_ids
	vector_unique( cut_t_ids );
	std::vector<int> tmp;
	std::set_difference( subdivide_t_ids.begin(), subdivide_t_ids.end(), cut_t_ids.begin(), cut_t_ids.end(), std::back_inserter( tmp ) );
	std::vector<bool> is_mark_surface( cut_t_ids.size(), true );
	cut_t_ids.insert( cut_t_ids.end(), tmp.begin(), tmp.end() );
	is_mark_surface.resize( is_mark_surface.size() + tmp.size(), false );

	std::vector<MeshTet> new_tets;
	std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
	std::vector<int> modified_t_ids;
	if( !subdivide_tets( insert_f_id, mesh, cut_mesh, points, map_edge_to_intersecting_point, track_surface_fs, cut_t_ids, is_mark_surface, new_tets,
		  new_track_surface_fs, modified_t_ids ) )
	{
		if( is_again )
		{
			if( is_uninserted_face_covered( insert_f_id, input_vertices, input_faces, cut_t_ids, mesh ) )
				return true;
		}
		mesh.logger().logWarning( "FAIL subdivide_tets" );
		return false;
	}

	push_new_tets( mesh, track_surface_fs, points, new_tets, new_track_surface_fs, modified_t_ids, is_again );

	simplify_subdivision_result( insert_f_id, input_vertices.size(), mesh, tree, track_surface_fs );
	return true;
}

void floatTetWild::push_new_tets( Mesh& mesh, std::vector<std::array<std::vector<int>, 4>>& track_surface_fs, std::vector<Vector3>& points,
  std::vector<MeshTet>& new_tets, std::vector<std::array<std::vector<int>, 4>>& new_track_surface_fs, std::vector<int>& modified_t_ids, bool is_again )
{
	const int old_v_size = mesh.tet_vertices.size();
	mesh.tet_vertices.resize( mesh.tet_vertices.size() + points.size() );
	for( int i = 0; i < points.size(); i++ )
		mesh.tet_vertices[ old_v_size + i ].pos = points[ i ];

	for( int i = 0; i < new_tets.size(); i++ )
	{
		if( i < modified_t_ids.size() )
		{
			for( int j = 0; j < 4; j++ )
				mesh.tet_vertices[ mesh.tets[ modified_t_ids[ i ] ][ j ] ].connTets.remove( modified_t_ids[ i ] );

			mesh.tets[ modified_t_ids[ i ] ] = new_tets[ i ];
			track_surface_fs[ modified_t_ids[ i ] ] = new_track_surface_fs[ i ];
			for( int j = 0; j < 4; j++ )
				mesh.tet_vertices[ mesh.tets[ modified_t_ids[ i ] ][ j ] ].connTets.add( modified_t_ids[ i ] );
		}
		else
		{
			for( int j = 0; j < 4; j++ )
				mesh.tet_vertices[ new_tets[ i ][ j ] ].connTets.add( mesh.tets.size() + i - modified_t_ids.size() );
		}
	}
	mesh.tets.insert( mesh.tets.end(), new_tets.begin() + modified_t_ids.size(), new_tets.end() );
	track_surface_fs.insert( track_surface_fs.end(), new_track_surface_fs.begin() + modified_t_ids.size(), new_track_surface_fs.end() );
	modified_t_ids.clear();
}

#include "EdgeCollapsing.h"

void floatTetWild::simplify_subdivision_result(
  int insert_f_id, int input_v_size, Mesh& mesh, AABBWrapper& tree, std::vector<std::array<std::vector<int>, 4>>& track_surface_fs )
{
	std::vector<std::array<int, 3>>& covered_tet_fs = mesh.globalVars.triangleInsertion.covered_tet_fs;
	if( covered_tet_fs.empty() )
		return;

	for( int i = 0; i < covered_tet_fs.size(); i++ )
		std::sort( covered_tet_fs[ i ].begin(), covered_tet_fs[ i ].end() );
	vector_unique( covered_tet_fs );

	std::vector<std::array<int, 3>> edges;
	for( int i = 0; i < covered_tet_fs.size(); i++ )
	{
		const auto f = covered_tet_fs[ i ];
		for( int j = 0; j < 3; j++ )
		{
			if( f[ j ] < f[ ( j + 1 ) % 3 ] )
				edges.push_back( { { f[ j ], f[ ( j + 1 ) % 3 ], i } } );
			else
				edges.push_back( { { f[ ( j + 1 ) % 3 ], f[ j ], i } } );
		}
	}
	std::sort( edges.begin(), edges.end(),
	  []( const std::array<int, 3>& a, const std::array<int, 3>& b ) { return std::make_tuple( a[ 0 ], a[ 1 ] ) < std::make_tuple( b[ 0 ], b[ 1 ] ); } );
	//
	std::unordered_set<int> freezed_v_ids;
	bool is_duplicated = false;
	for( int i = 0; i < edges.size() - 1; i++ )
	{
		if( edges[ i ][ 0 ] < input_v_size )
			freezed_v_ids.insert( edges[ i ][ 0 ] );
		if( edges[ i ][ 1 ] < input_v_size )
			freezed_v_ids.insert( edges[ i ][ 1 ] );

		if( edges[ i ][ 0 ] == edges[ i + 1 ][ 0 ] && edges[ i ][ 1 ] == edges[ i + 1 ][ 1 ] )
		{
			is_duplicated = true;
			edges.erase( edges.begin() + i );
			i--;
		}
		else
		{
			if( !is_duplicated )
			{
				// boundary edges
				freezed_v_ids.insert( edges[ i ][ 0 ] );
				freezed_v_ids.insert( edges[ i ][ 1 ] );
			}
			is_duplicated = false;
		}
	}
	if( !is_duplicated )
	{
		freezed_v_ids.insert( edges.back()[ 0 ] );
		freezed_v_ids.insert( edges.back()[ 1 ] );
	}

	std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_s> ec_queue;
	for( const auto& e : edges )
	{
		Scalar l_2 = get_edge_length_2( mesh, e[ 0 ], e[ 1 ] );
		if( freezed_v_ids.find( e[ 0 ] ) == freezed_v_ids.end() )
			ec_queue.push( ElementInQueue( { { e[ 0 ], e[ 1 ] } }, l_2 ) );
		if( freezed_v_ids.find( e[ 1 ] ) == freezed_v_ids.end() )
			ec_queue.push( ElementInQueue( { { e[ 1 ], e[ 0 ] } }, l_2 ) );
	}
	if( ec_queue.empty() )
		return;

	std::unordered_set<int> all_v_ids;
	for( const auto& e : edges )
	{
		all_v_ids.insert( e[ 0 ] );
		all_v_ids.insert( e[ 1 ] );
	}

	int _ts = 0;
	std::vector<int> _tet_tss;
	bool is_update_tss = false;
	int cnt_suc = 0;
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

		if( !is_valid_edge( mesh, v_ids[ 0 ], v_ids[ 1 ] ) )
			continue;
		if( freezed_v_ids.find( v_ids[ 0 ] ) != freezed_v_ids.end() )
			continue;
		Scalar weight = get_edge_length_2( mesh, v_ids[ 0 ], v_ids[ 1 ] );
		if( weight != old_weight )
			continue;

		// check track_surface_fs
		int v1_id = v_ids[ 0 ];
		int v2_id = v_ids[ 1 ];
		bool is_valid = true;
		for( int t_id : mesh.tet_vertices[ v1_id ].connTets )
		{
			for( int j = 0; j < 4; j++ )
			{
				if( ( !track_surface_fs[ t_id ][ j ].empty() && !vector_contains( track_surface_fs[ t_id ][ j ], insert_f_id ) ) ||
					( mesh.tets[ t_id ][ j ] != v1_id &&
					  ( mesh.tets[ t_id ].is_surface_fs[ j ] != NOT_SURFACE || mesh.tets[ t_id ].is_bbox_fs[ j ] != NOT_BBOX ) ) )
				{
					is_valid = false;
					break;
				}
			}
			if( !is_valid )
				break;
		}
		if( !is_valid )
			continue;

		EdgesSet new_edges;
		static const bool is_check_quality = true;
		const auto& v1_conn_tets = mesh.tet_vertices[ v1_id ].connTets;
		for( int t_id : v1_conn_tets )
			mesh.tets[ t_id ].quality = get_quality( mesh, t_id );

		const eCollapseStatus result = collapse_an_edge( mesh, v_ids[ 0 ], v_ids[ 1 ], tree, new_edges, _ts, _tet_tss, is_check_quality, is_update_tss );
		if( isSuccessStatus( result ) )
		{
			new_edges.enumerate(
			  [ & ]( int e0, int e1 )
			  {
				  if( all_v_ids.find( e0 ) == all_v_ids.end() || all_v_ids.find( e1 ) == all_v_ids.end() )
					  return;
				  Scalar l_2 = get_edge_length_2( mesh, e0, e1 );
				  if( freezed_v_ids.find( e0 ) == freezed_v_ids.end() )
					  ec_queue.push( ElementInQueue( e0, e1, l_2 ) );
				  if( freezed_v_ids.find( e1 ) == freezed_v_ids.end() )
					  ec_queue.push( ElementInQueue( e1, e0, l_2 ) );
			  } );
			cnt_suc++;
		}
	}
}

namespace
{
	using namespace floatTetWild;

	static void findCuttingTetsOmp( Mesh& mesh, __m256d min_f, __m256d max_f, std::queue<int>& queue_t_ids, std::vector<bool>& is_visited )
	{
		const int length = (int)mesh.tets.size();
#pragma omp parallel for
		for( int t_id = 0; t_id < mesh.tets.size(); t_id++ )
		{
			const auto& mt = mesh.tets[ t_id ];
			if( mt.is_removed )
				continue;

			__m256d min_t, max_t;
			get_bbox_tet( mesh.tet_vertices[ mt[ 0 ] ].pos, mesh.tet_vertices[ mt[ 1 ] ].pos, mesh.tet_vertices[ mt[ 2 ] ].pos,
			  mesh.tet_vertices[ mt[ 3 ] ].pos, min_t, max_t );

			if( !is_bbox_intersected( min_f, max_f, min_t, max_t ) )
				continue;

#pragma omp critical
			{
				queue_t_ids.push( t_id );
				is_visited[ t_id ] = true;
			}
		}
	}

	static void findCuttingTets( Mesh& mesh, __m256d min_f, __m256d max_f, std::queue<int>& queue_t_ids, std::vector<bool>& is_visited )
	{
		for( size_t t_id = 0; t_id < mesh.tets.size(); t_id++ )
		{
			const auto& mt = mesh.tets[ t_id ];
			if( mt.is_removed )
				continue;

			__m256d min_t, max_t;
			get_bbox_tet( mesh.tet_vertices[ mt[ 0 ] ].pos, mesh.tet_vertices[ mt[ 1 ] ].pos, mesh.tet_vertices[ mt[ 2 ] ].pos,
			  mesh.tet_vertices[ mt[ 3 ] ].pos, min_t, max_t );

			if( !is_bbox_intersected( min_f, max_f, min_t, max_t ) )
				continue;

			queue_t_ids.push( t_id );
			is_visited[ t_id ] = true;
		}
	}
}  // namespace

void floatTetWild::find_cutting_tets( int f_id, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  const std::array<Vector3, 3>& vs, Mesh& mesh, std::vector<int>& cut_t_ids, bool is_again )
{
	auto tm = mesh.times.findCuttingTets.measure();

	FindCuttingTetsBuffers& buffers = mesh.findCuttingTetsBuffers;

	std::vector<bool>& is_visited = buffers.is_visited;
	is_visited.clear();
	is_visited.resize( mesh.tets.size(), false );

	std::queue<int>& queue_t_ids = buffers.queue_t_ids;

	if( !is_again )
	{
		std::vector<int>& n_t_ids = buffers.n_t_ids;
		n_t_ids.clear();

		for( int j = 0; j < 3; j++ )
		{
			const auto& conn_tets = mesh.tet_vertices[ input_faces[ f_id ][ j ] ].connTets;
			n_t_ids.insert( n_t_ids.end(), conn_tets.begin(), conn_tets.end() );
		}
		vector_unique( n_t_ids );

		for( int t_id : n_t_ids )
		{
			is_visited[ t_id ] = true;
			queue_t_ids.push( t_id );
		}
	}
	else
	{
		__m256d min_f, max_f;
		get_bbox_face(
		  input_vertices[ input_faces[ f_id ][ 0 ] ], input_vertices[ input_faces[ f_id ][ 1 ] ], input_vertices[ input_faces[ f_id ][ 2 ] ], min_f, max_f );
		if( mesh.params.num_threads > 1 )
			findCuttingTetsOmp( mesh, min_f, max_f, queue_t_ids, is_visited );
		else
			findCuttingTets( mesh, min_f, max_f, queue_t_ids, is_visited );
	}

	while( !queue_t_ids.empty() )
	{
		int t_id = queue_t_ids.front();
		queue_t_ids.pop();

		if( is_again )
		{
			bool is_cut = false;
			for( int j = 0; j < 3; j++ )
			{
				if( is_point_inside_tet( input_vertices[ input_faces[ f_id ][ j ] ], mesh.tet_vertices[ mesh.tets[ t_id ][ 0 ] ].pos,
					  mesh.tet_vertices[ mesh.tets[ t_id ][ 1 ] ].pos, mesh.tet_vertices[ mesh.tets[ t_id ][ 2 ] ].pos,
					  mesh.tet_vertices[ mesh.tets[ t_id ][ 3 ] ].pos ) )
				{
					is_cut = true;
					break;
				}
			}
			if( is_cut )
			{
				cut_t_ids.push_back( t_id );
				for( int j = 0; j < 4; j++ )
				{
					for( int n_t_id : mesh.tet_vertices[ mesh.tets[ t_id ][ j ] ].connTets )
					{
						if( !is_visited[ n_t_id ] )
						{
							is_visited[ n_t_id ] = true;
							queue_t_ids.push( n_t_id );
						}
					}
				}
				continue;
			}
		}

		std::array<eOrientation, 4> oris;
		for( int j = 0; j < 4; j++ )
			oris[ j ] = Predicates::orient_3d( vs[ 0 ], vs[ 1 ], vs[ 2 ], mesh.tet_vertices[ mesh.tets[ t_id ][ j ] ].pos );

		bool is_cutted = false;
		std::array<bool, 4> is_cut_vs = { { false, false, false, false } };	 /// is v on cut face
		for( int j = 0; j < 4; j++ )
		{
			int cnt_pos = 0;
			int cnt_neg = 0;
			int cnt_on = 0;
			for( int k = 0; k < 3; k++ )
			{
				if( oris[ ( j + k + 1 ) % 4 ] == eOrientation::Zero )
					cnt_on++;
				else if( oris[ ( j + k + 1 ) % 4 ] == eOrientation::Positive )
					cnt_pos++;
				else
					cnt_neg++;
			}

			eCutResult result = eCutResult::Empty;
			auto& tp1 = mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 1 ) % 4 ] ].pos;
			auto& tp2 = mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 2 ) % 4 ] ].pos;
			auto& tp3 = mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 3 ) % 4 ] ].pos;

			if( cnt_on == 3 )
			{
				if( is_tri_tri_cutted_hint( vs[ 0 ], vs[ 1 ], vs[ 2 ], tp1, tp2, tp3, eCutResult::Coplanar ) == eCutResult::Coplanar )
				{
					result = eCutResult::Coplanar;
					is_cutted = true;
					is_cut_vs[ ( j + 1 ) % 4 ] = true;
					is_cut_vs[ ( j + 2 ) % 4 ] = true;
					is_cut_vs[ ( j + 3 ) % 4 ] = true;
				}
			}
			else if( cnt_pos > 0 && cnt_neg > 0 )
			{
				if( is_tri_tri_cutted_hint( vs[ 0 ], vs[ 1 ], vs[ 2 ], tp1, tp2, tp3, eCutResult::Face ) == eCutResult::Face )
				{
					result = eCutResult::Face;
					is_cutted = true;
					is_cut_vs[ ( j + 1 ) % 4 ] = true;
					is_cut_vs[ ( j + 2 ) % 4 ] = true;
					is_cut_vs[ ( j + 3 ) % 4 ] = true;
				}
			}
			else if( cnt_on == 2 && oris[ ( j + 1 ) % 4 ] == eOrientation::Zero && oris[ ( j + 2 ) % 4 ] == eOrientation::Zero )
			{
				if( is_tri_tri_cutted_hint( vs[ 0 ], vs[ 1 ], vs[ 2 ], tp1, tp2, tp3, eCutResult::Edge0 ) == eCutResult::Edge0 )
				{
					result = eCutResult::Edge0;
					is_cut_vs[ ( j + 1 ) % 4 ] = true;
					is_cut_vs[ ( j + 2 ) % 4 ] = true;
				}
			}
			else if( cnt_on == 2 && oris[ ( j + 2 ) % 4 ] == eOrientation::Zero && oris[ ( j + 3 ) % 4 ] == eOrientation::Zero )
			{
				if( is_tri_tri_cutted_hint( vs[ 0 ], vs[ 1 ], vs[ 2 ], tp1, tp2, tp3, eCutResult::Edge1 ) == eCutResult::Edge1 )
				{
					result = eCutResult::Edge1;
					is_cut_vs[ ( j + 2 ) % 4 ] = true;
					is_cut_vs[ ( j + 3 ) % 4 ] = true;
				}
			}
			else if( cnt_on == 2 && oris[ ( j + 3 ) % 4 ] == eOrientation::Zero && oris[ ( j + 1 ) % 4 ] == eOrientation::Zero )
			{
				if( is_tri_tri_cutted_hint( vs[ 0 ], vs[ 1 ], vs[ 2 ], tp1, tp2, tp3, eCutResult::Edge2 ) == eCutResult::Edge2 )
				{
					result = eCutResult::Edge2;
					is_cut_vs[ ( j + 3 ) % 4 ] = true;
					is_cut_vs[ ( j + 1 ) % 4 ] = true;
				}
			}
		}
		if( is_cutted )
			cut_t_ids.push_back( t_id );

		for( int j = 0; j < 4; j++ )
		{
			if( !is_cut_vs[ j ] )
				continue;
			for( int n_t_id : mesh.tet_vertices[ mesh.tets[ t_id ][ j ] ].connTets )
			{
				if( !is_visited[ n_t_id ] )
				{
					is_visited[ n_t_id ] = true;
					queue_t_ids.push( n_t_id );
				}
			}
		}
	}
}

namespace
{
	static const std::array<std::array<int, 2>, 6> t_es = { { { { 0, 1 } }, { { 1, 2 } }, { { 2, 0 } }, { { 0, 3 } }, { { 1, 3 } }, { { 2, 3 } } } };
	static const std::array<std::array<int, 3>, 4> t_f_es = { { { { 1, 5, 4 } }, { { 5, 3, 2 } }, { { 3, 0, 4 } }, { { 0, 1, 2 } } } };
	static const std::array<std::array<int, 3>, 4> t_f_vs = { { { { 3, 1, 2 } }, { { 0, 2, 3 } }, { { 1, 3, 0 } }, { { 2, 0, 1 } } } };
}

bool floatTetWild::subdivide_tets( int insert_f_id, Mesh& mesh, CutMesh& cut_mesh, std::vector<Vector3>& points,
  std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point, std::vector<std::array<std::vector<int>, 4>>& track_surface_fs,
  std::vector<int>& subdivide_t_ids, std::vector<bool>& is_mark_surface, std::vector<MeshTet>& new_tets,
  std::vector<std::array<std::vector<int>, 4>>& new_track_surface_fs, std::vector<int>& modified_t_ids )
{
	auto tm = mesh.times.subdivideTets.measure();

	std::vector<std::array<int, 3>>& covered_tet_fs = mesh.globalVars.triangleInsertion.covered_tet_fs;
	covered_tet_fs.clear();

	for( int I = 0; I < subdivide_t_ids.size(); I++ )
	{
		int t_id = subdivide_t_ids[ I ];
		bool is_mark_sf = is_mark_surface[ I ];

		/////
		std::bitset<6> config_bit;
		std::array<std::pair<int, int>, 6> on_edge_p_ids;
		int cnt = 4;
		for( int i = 0; i < t_es.size(); i++ )
		{
			const auto& le = t_es[ i ];
			std::array<int, 2> e = { { mesh.tets[ t_id ][ le[ 0 ] ], mesh.tets[ t_id ][ le[ 1 ] ] } };
			sortInt2( e );
			if( map_edge_to_intersecting_point.find( e ) == map_edge_to_intersecting_point.end() )
			{
				on_edge_p_ids[ i ].first = -1;
				on_edge_p_ids[ i ].second = -1;
			}
			else
			{
				on_edge_p_ids[ i ].first = cnt++;
				on_edge_p_ids[ i ].second = map_edge_to_intersecting_point[ e ];
				config_bit.set( i );
			}
		}
		int config_id = config_bit.to_ulong();
		if( config_id == 0 )
		{  // no intersection
			if( is_mark_sf )
			{
				for( int j = 0; j < 4; j++ )
				{
					int cnt_on = 0;
					for( int k = 0; k < 3; k++ )
					{
						myassert( cut_mesh.map_v_ids.find( mesh.tets[ t_id ][ ( j + k + 1 ) % 4 ] ) != cut_mesh.map_v_ids.end(),
						  "cut_mesh.map_v_ids.find(mesh.tets[t_id][(j + k + 1) % 4]) != cut_mesh.map_v_ids.end()!!" );	// fortest
						if( cut_mesh.is_v_on_plane( cut_mesh.map_v_ids[ mesh.tets[ t_id ][ ( j + k + 1 ) % 4 ] ] ) )
						{
							cnt_on++;
						}
					}
					// fortest
					if( cnt_on == 4 )
						mesh.logger().logWarning( "cnt_on==4!!" );
					// fortest

					if( cnt_on == 3 )
					{
						new_tets.push_back( mesh.tets[ t_id ] );
						new_track_surface_fs.push_back( track_surface_fs[ t_id ] );
						( new_track_surface_fs.back() )[ j ].push_back( insert_f_id );
						modified_t_ids.push_back( t_id );

						covered_tet_fs.push_back(
						  { { mesh.tets[ t_id ][ ( j + 1 ) % 4 ], mesh.tets[ t_id ][ ( j + 2 ) % 4 ], mesh.tets[ t_id ][ ( j + 3 ) % 4 ] } } );
						break;
					}
				}
			}
			continue;
		}

		const auto& configs = CutTable::get_tet_confs( config_id );
		if( configs.empty() )
			continue;

		/////
		SmallBuffer<Vector2i, 4> my_diags;
		for( int j = 0; j < 4; j++ )
		{
			SmallBuffer<int, 3> le_ids;
			for( int k = 0; k < 3; k++ )
			{
				if( on_edge_p_ids[ t_f_es[ j ][ k ] ].first < 0 )
					continue;
				le_ids.push_back( k );
			}
			if( le_ids.size() != 2 )  // no ambiguity
				continue;

			auto& diag = my_diags.emplace_back();
			if( on_edge_p_ids[ t_f_es[ j ][ le_ids[ 0 ] ] ].second > on_edge_p_ids[ t_f_es[ j ][ le_ids[ 1 ] ] ].second )
				diag << on_edge_p_ids[ t_f_es[ j ][ le_ids[ 0 ] ] ].first, t_f_vs[ j ][ le_ids[ 0 ] ];
			else
				diag << on_edge_p_ids[ t_f_es[ j ][ le_ids[ 1 ] ] ].first, t_f_vs[ j ][ le_ids[ 1 ] ];

			sortInt2( diag );
		}
		std::sort( my_diags.begin(), my_diags.end(), []( const Vector2i& a, const Vector2i& b ) { return compareInt2( a, b ); } );

		/////
		std::map<int, int> map_lv_to_v_id;
		const int v_size = mesh.tet_vertices.size();
		const int vp_size = mesh.tet_vertices.size() + points.size();
		for( int i = 0; i < 4; i++ )
			map_lv_to_v_id[ i ] = mesh.tets[ t_id ][ i ];
		cnt = 0;
		for( int i = 0; i < t_es.size(); i++ )
		{
			if( config_bit[ i ] == 0 )
				continue;
			map_lv_to_v_id[ 4 + cnt ] = v_size + on_edge_p_ids[ i ].second;
			cnt++;
		}

		/////
		auto get_centroid = [ & ]( const CutTable::Vec4Buffer& config, int lv_id, Vector3& c )
		{
			std::vector<int> n_ids;
			for( const auto& tet : config )
			{
				std::vector<int> tmp;
				for( int j = 0; j < 4; j++ )
				{
					if( tet[ j ] != lv_id )
						tmp.push_back( tet[ j ] );
				}
				if( tmp.size() == 4 )
					continue;
				n_ids.insert( n_ids.end(), tmp.begin(), tmp.end() );
			}
			vector_unique( n_ids );
			c << 0, 0, 0;
			for( int n_id : n_ids )
			{
				int v_id = map_lv_to_v_id[ n_id ];
				if( v_id < v_size )
					c += mesh.tet_vertices[ v_id ].pos;
				else
					c += points[ v_id - v_size ];
			}
			c /= n_ids.size();
		};

		auto check_config = [ & ]( int diag_config_id, std::vector<std::pair<int, Vector3>>& centroids )
		{
			const auto& config = CutTable::get_tet_conf( config_id, diag_config_id );
			Scalar min_q = -666;
			int cnt = 0;
			std::map<int, int> map_lv_to_c;
			for( const auto& tet : config )
			{
				std::array<Vector3, 4> vs;
				for( int j = 0; j < 4; j++ )
				{
					if( map_lv_to_v_id.find( tet[ j ] ) == map_lv_to_v_id.end() )
					{
						if( map_lv_to_c.find( tet[ j ] ) == map_lv_to_c.end() )
						{
							get_centroid( config, tet[ j ], vs[ j ] );
							map_lv_to_c[ tet[ j ] ] = centroids.size();
							centroids.push_back( std::make_pair( tet[ j ], vs[ j ] ) );
						}
						vs[ j ] = centroids[ map_lv_to_c[ tet[ j ] ] ].second;
					}
					else
					{
						int v_id = map_lv_to_v_id[ tet[ j ] ];
						if( v_id < v_size )
							vs[ j ] = mesh.tet_vertices[ v_id ].pos;
						else
							vs[ j ] = points[ v_id - v_size ];
					}
				}

				Scalar volume = Predicates::orient_3d_volume( vs[ 0 ], vs[ 1 ], vs[ 2 ], vs[ 3 ] );

				//                //fortest
				//                if(volume==0) {
				//                    cout<<std::setprecision(16)<<endl;
				//                    cout << "volume = " << volume << ",ori = " << Predicates::orient_3d(vs[0], vs[1], vs[2], vs[3]) << endl;
				//                    cout << "centroids.size = " << centroids.size() << endl;
				//                    cout << "config_id = " << config_id << endl;
				//
				//                    cout<<"vertices"<<endl;
				//                    for(int k=0;k<4;k++){
				//                        cout<<"v"<<mesh.tets[t_id][k]<<": "<<mesh.tet_vertices[mesh.tets[t_id][k]].pos.transpose()<<endl;
				//                    }
				//                    cout<<"intersecting points"<<endl;
				//                    for (int i = 0; i < t_es.size(); i++) {
				//                        const auto &le = t_es[i];
				//                        std::array<int, 2> e = {{mesh.tets[t_id][le[0]], mesh.tets[t_id][le[1]]}};
				//                        if (e[0] > e[1])
				//                            std::swap(e[0], e[1]);
				//                        if (map_edge_to_intersecting_point.find(e) != map_edge_to_intersecting_point.end()) {
				//                            cout << e[0] << " " << e[1] << ": p" << map_edge_to_intersecting_point[e] << " "
				//                                 << points[map_edge_to_intersecting_point[e]].transpose() << endl;
				//                            cout<<"is_mark_sf = "<<is_mark_sf<<endl;
				//                            if(cut_mesh.map_v_ids.find(e[0])!=cut_mesh.map_v_ids.end()) {
				//                                cout << "e[0] to_plane_dists = " << cut_mesh.to_plane_dists[cut_mesh.map_v_ids[e[0]]]
				//                                     << endl;
				//                                cout << "e[0] is_snapped = " << cut_mesh.is_snapped[cut_mesh.map_v_ids[e[0]]]
				//                                     << endl;
				//                            }
				//                            if(cut_mesh.map_v_ids.find(e[1])!=cut_mesh.map_v_ids.end()) {
				//                                cout << "e[1] to_plane_dists = " << cut_mesh.to_plane_dists[cut_mesh.map_v_ids[e[1]]]
				//                                     << endl;
				//                                cout << "e[1] is_snapped = " << cut_mesh.is_snapped[cut_mesh.map_v_ids[e[1]]]
				//                                     << endl;
				//                            }
				//
				////                            Vector3 p;
				////                            Scalar _;
				////                            seg_plane_intersection(mesh.tet_vertices[e[0]].pos, mesh.tet_vertices[e[1]].pos,
				////                                                   cut_mesh.p_vs[0], cut_mesh.p_n, p, _);
				////                            cout << "new p:" << p.transpose() << endl;
				////                            cout << (p - mesh.tet_vertices[e[1]].pos).squaredNorm() << " " << SCALAR_ZERO_2 << endl;
				//                        }
				//                    }
				//
				//                    pausee();
				//                }
				//                //fortest

				if( cnt == 0 )
					min_q = volume;
				else if( volume < min_q )
					min_q = volume;
				cnt++;
			}

			return min_q;
		};

		int diag_config_id = 0;
		std::vector<std::pair<int, Vector3>> centroids;
		if( !my_diags.empty() )
		{
			const auto& all_diags = CutTable::get_diag_confs( config_id );
			std::vector<std::pair<int, Scalar>> min_qualities( all_diags.size() );
			std::vector<std::vector<std::pair<int, Vector3>>> all_centroids( all_diags.size() );
			for( int i = 0; i < all_diags.size(); i++ )
			{
				if( !all_diags[ i ].equal( my_diags.data(), my_diags.data() + my_diags.size() ) )
				{
					min_qualities[ i ] = std::make_pair( i, -1 );
					continue;
				}

				std::vector<std::pair<int, Vector3>> tmp_centroids;
				Scalar min_q = check_config( i, tmp_centroids );
				//                if (min_q < SCALAR_ZERO_3)
				//                    continue;
				min_qualities[ i ] = std::make_pair( i, min_q );
				all_centroids[ i ] = tmp_centroids;
			}
			std::sort( min_qualities.begin(), min_qualities.end(),
			  []( const std::pair<int, Scalar>& a, const std::pair<int, Scalar>& b ) { return a.second < b.second; } );

			if( min_qualities.back().second < SCALAR_ZERO_3 )
			{  // if tet quality is too bad
				//                cout<<std::setprecision(16)<<"return 1 "<<min_qualities.back().second<<endl;
				return false;
			}

			diag_config_id = min_qualities.back().first;
			centroids = all_centroids[ diag_config_id ];
		}
		else
		{
			Scalar min_q = check_config( diag_config_id, centroids );
			if( min_q < SCALAR_ZERO_3 )
			{
				//                cout<<std::setprecision(16)<<"return 2 "<<min_q<<endl;
				return false;
			}
		}

		for( int i = 0; i < centroids.size(); i++ )
		{
			map_lv_to_v_id[ centroids[ i ].first ] = vp_size + i;
			points.push_back( centroids[ i ].second );
		}

		// add new tets
		const auto& config = CutTable::get_tet_conf( config_id, diag_config_id );
		const auto& new_is_surface_fs = CutTable::get_surface_conf( config_id, diag_config_id );
		const auto& new_local_f_ids = CutTable::get_face_id_conf( config_id, diag_config_id );
		for( int i = 0; i < config.size(); i++ )
		{
			const auto& t = config[ i ];
			new_tets.push_back( MeshTet( map_lv_to_v_id[ t[ 0 ] ], map_lv_to_v_id[ t[ 1 ] ], map_lv_to_v_id[ t[ 2 ] ], map_lv_to_v_id[ t[ 3 ] ] ) );

			// fortest
			//            cout << map_lv_to_v_id[t[0]] << " " << map_lv_to_v_id[t[1]] << " " << map_lv_to_v_id[t[2]] << " "
			//                 << map_lv_to_v_id[t[3]] << endl;

			new_track_surface_fs.emplace_back();
			for( int j = 0; j < 4; j++ )
			{
				if( new_is_surface_fs[ i ][ j ] && is_mark_sf )
				{
					( new_track_surface_fs.back() )[ j ].push_back( insert_f_id );

					covered_tet_fs.push_back( { { new_tets.back()[ ( j + 1 ) % 4 ], new_tets.back()[ ( j + 3 ) % 4 ], new_tets.back()[ ( j + 2 ) % 4 ] } } );
				}

				int old_local_f_id = new_local_f_ids[ i ][ j ];
				if( old_local_f_id < 0 )
					continue;
				( new_track_surface_fs.back() )[ j ].insert( ( new_track_surface_fs.back() )[ j ].end(), track_surface_fs[ t_id ][ old_local_f_id ].begin(),
				  track_surface_fs[ t_id ][ old_local_f_id ].end() );
				( new_tets.back() ).is_bbox_fs[ j ] = mesh.tets[ t_id ].is_bbox_fs[ old_local_f_id ];
				( new_tets.back() ).is_surface_fs[ j ] = mesh.tets[ t_id ].is_surface_fs[ old_local_f_id ];
				( new_tets.back() ).surface_tags[ j ] = mesh.tets[ t_id ].surface_tags[ old_local_f_id ];
			}
		}
		modified_t_ids.push_back( t_id );
	}

	return true;
}

void floatTetWild::pair_track_surface_fs( Mesh& mesh, std::vector<std::array<std::vector<int>, 4>>& track_surface_fs )
{
	std::vector<std::array<bool, 4>> is_visited( track_surface_fs.size(), { { false, false, false, false } } );
	for( int t_id = 0; t_id < track_surface_fs.size(); t_id++ )
	{
		for( int j = 0; j < 4; j++ )
		{
			//
			if( is_visited[ t_id ][ j ] )
				continue;
			is_visited[ t_id ][ j ] = true;
			if( track_surface_fs[ t_id ][ j ].empty() )
				continue;
			//
			int opp_t_id = get_opp_t_id( t_id, j, mesh );
			if( opp_t_id < 0 )
				continue;
			int k =
			  get_local_f_id( opp_t_id, mesh.tets[ t_id ][ ( j + 1 ) % 4 ], mesh.tets[ t_id ][ ( j + 2 ) % 4 ], mesh.tets[ t_id ][ ( j + 3 ) % 4 ], mesh );
			is_visited[ opp_t_id ][ k ] = true;
			//
			std::sort( track_surface_fs[ t_id ][ j ].begin(), track_surface_fs[ t_id ][ j ].end() );
			std::sort( track_surface_fs[ opp_t_id ][ k ].begin(), track_surface_fs[ opp_t_id ][ k ].end() );
			std::vector<int> f_ids;
			if( track_surface_fs[ t_id ][ j ] != track_surface_fs[ opp_t_id ][ k ] )
			{
				std::set_union( track_surface_fs[ t_id ][ j ].begin(), track_surface_fs[ t_id ][ j ].end(), track_surface_fs[ opp_t_id ][ k ].begin(),
				  track_surface_fs[ opp_t_id ][ k ].end(), std::back_inserter( f_ids ) );
				track_surface_fs[ t_id ][ j ] = f_ids;
				track_surface_fs[ opp_t_id ][ k ] = f_ids;
			}
		}
	}
}

void floatTetWild::find_boundary_edges( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  const std::vector<bool>& is_face_inserted, const std::vector<bool>& old_is_face_inserted,
  std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos, std::vector<bool>& is_on_cut_edges, std::vector<std::array<int, 2>>& b_edges,
  const Logger& log )
{
	std::vector<std::array<int, 2>> edges;
	std::vector<std::vector<int>> conn_tris( input_vertices.size() );
	std::vector<std::vector<int>> uninserted_conn_tris( input_vertices.size() );
	for( int i = 0; i < input_faces.size(); i++ )
	{
		if( !is_face_inserted[ i ] )
		{  /// use currently inserted faces as mesh
			for( int j = 0; j < 3; j++ )
				uninserted_conn_tris[ input_faces[ i ][ j ] ].push_back( i );
			continue;
		}
		const auto& f = input_faces[ i ];
		for( int j = 0; j < 3; j++ )
		{
			// edges
			std::array<int, 2> e = { { f[ j ], f[ ( j + 1 ) % 3 ] } };
			if( e[ 0 ] > e[ 1 ] )
				std::swap( e[ 0 ], e[ 1 ] );
			edges.push_back( e );
			// conn_tris
			conn_tris[ input_faces[ i ][ j ] ].push_back( i );
		}
	}
	vector_unique( edges );

	int cnt1 = 0;
	int cnt2 = 0;
	for( const auto& e : edges )
	{
		std::vector<int> n12_f_ids;
		std::set_intersection(
		  conn_tris[ e[ 0 ] ].begin(), conn_tris[ e[ 0 ] ].end(), conn_tris[ e[ 1 ] ].begin(), conn_tris[ e[ 1 ] ].end(), std::back_inserter( n12_f_ids ) );
		std::vector<int> uninserted_n12_f_ids;
		std::set_intersection( uninserted_conn_tris[ e[ 0 ] ].begin(), uninserted_conn_tris[ e[ 0 ] ].end(), uninserted_conn_tris[ e[ 1 ] ].begin(),
		  uninserted_conn_tris[ e[ 1 ] ].end(), std::back_inserter( uninserted_n12_f_ids ) );

		bool needs_preserve = false;
		for( int f_id : n12_f_ids )
		{
			if( !old_is_face_inserted[ f_id ] )
			{
				needs_preserve = true;
				break;
			}
		}

		if( n12_f_ids.size() == 1 )
		{  // open boundary
			b_edges.push_back( e );
			if( needs_preserve )
			{
				b_edge_infos.push_back( std::make_pair( e, n12_f_ids ) );
				if( !uninserted_n12_f_ids.empty() )
					is_on_cut_edges.push_back( true );
				else
					is_on_cut_edges.push_back( false );
			}
			cnt1++;
		}
		else
		{
			int f_id = n12_f_ids[ 0 ];
			int j = 0;
			for( ; j < 3; j++ )
			{
				if( ( input_faces[ f_id ][ j ] == e[ 0 ] && input_faces[ f_id ][ mod3( j + 1 ) ] == e[ 1 ] ) ||
					( input_faces[ f_id ][ j ] == e[ 1 ] && input_faces[ f_id ][ mod3( j + 1 ) ] == e[ 0 ] ) )
					break;
			}

			Vector3 n =
			  get_normal( input_vertices[ input_faces[ f_id ][ 0 ] ], input_vertices[ input_faces[ f_id ][ 1 ] ], input_vertices[ input_faces[ f_id ][ 2 ] ] );
			int t = get_t( input_vertices[ input_faces[ f_id ][ 0 ] ], input_vertices[ input_faces[ f_id ][ 1 ] ], input_vertices[ input_faces[ f_id ][ 2 ] ] );

			bool is_fine = false;
			for( int k = 0; k < n12_f_ids.size(); k++ )
			{
				if( n12_f_ids[ k ] == f_id )
					continue;
				Vector3 n1 = get_normal( input_vertices[ input_faces[ n12_f_ids[ k ] ][ 0 ] ], input_vertices[ input_faces[ n12_f_ids[ k ] ][ 1 ] ],
				  input_vertices[ input_faces[ n12_f_ids[ k ] ][ 2 ] ] );
				if( abs( n1.dot( n ) ) < 1 - SCALAR_ZERO )
				{
					is_fine = true;
					break;
				}
			}
			if( is_fine )
				continue;

			is_fine = false;
			eOrientation ori;
			for( int k = 0; k < n12_f_ids.size(); k++ )
			{
				for( int r = 0; r < 3; r++ )
				{
					if( input_faces[ n12_f_ids[ k ] ][ r ] != input_faces[ f_id ][ j ] &&
						input_faces[ n12_f_ids[ k ] ][ r ] != input_faces[ f_id ][ mod3( j + 1 ) ] )
					{
						if( k == 0 )
						{
							ori = Predicates::orient_2d( to_2d( input_vertices[ input_faces[ f_id ][ j ] ], t ),
							  to_2d( input_vertices[ input_faces[ f_id ][ mod3( j + 1 ) ] ], t ),
							  to_2d( input_vertices[ input_faces[ n12_f_ids[ k ] ][ r ] ], t ) );
							break;
						}
						eOrientation new_ori = Predicates::orient_2d( to_2d( input_vertices[ input_faces[ f_id ][ j ] ], t ),
						  to_2d( input_vertices[ input_faces[ f_id ][ mod3( j + 1 ) ] ], t ),
						  to_2d( input_vertices[ input_faces[ n12_f_ids[ k ] ][ r ] ], t ) );
						if( new_ori != ori )
							is_fine = true;
						break;
					}
				}
				if( is_fine )
					break;
			}
			if( is_fine )
				continue;

			cnt2++;
			b_edges.push_back( e );
			if( needs_preserve )
			{
				b_edge_infos.push_back( std::make_pair( e, n12_f_ids ) );
				if( !uninserted_n12_f_ids.empty() )
					is_on_cut_edges.push_back( true );
				else
					is_on_cut_edges.push_back( false );
			}
		}
	}

	log.logDebug( "#boundary_e1 = %i, #boundary_e2 = %i", cnt1, cnt2 );
}

bool floatTetWild::insert_boundary_edges( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos, std::vector<bool>& is_on_cut_edges,
  std::vector<std::array<std::vector<int>, 4>>& track_surface_fs, Mesh& mesh, AABBWrapper& tree, std::vector<bool>& is_face_inserted, bool is_again,
  std::vector<std::array<int, 3>>& known_surface_fs, std::vector<std::array<int, 3>>& known_not_surface_fs )
{
	igl::Timer timer;

	auto mark_known_surface_fs = [ & ]( const std::array<int, 3>& f, int tag )
	{
		std::vector<int> n_t_ids;
		setIntersection( mesh.tet_vertices[ f[ 0 ] ].connTets, mesh.tet_vertices[ f[ 1 ] ].connTets, mesh.tet_vertices[ f[ 2 ] ].connTets, n_t_ids );
		if( n_t_ids.size() != 2 )  // todo:?????
			return;

		for( int t_id : n_t_ids )
		{
			int j = get_local_f_id( t_id, f[ 0 ], f[ 1 ], f[ 2 ], mesh );
			mesh.tets[ t_id ].is_surface_fs[ j ] = tag;
		}
	};

	auto record_boundary_info = [ & ]( const std::vector<Vector3>& points, const std::vector<int>& snapped_v_ids, const std::array<int, 2>& e, bool is_on_cut )
	{
		const int tet_vertices_size = mesh.tet_vertices.size();
		for( int i = points.size(); i > 0; i-- )
		{
			mesh.tet_vertices[ tet_vertices_size - i ].setFlag( eVertexFlags::Boundary );
			mesh.tet_vertices[ tet_vertices_size - i ].setFlag( eVertexFlags::Cut, is_on_cut );
		}

		for( int v_id : snapped_v_ids )
		{
			Scalar dis_2 = p_seg_squared_dist_3d( mesh.tet_vertices[ v_id ].pos, input_vertices[ e[ 0 ] ], input_vertices[ e[ 1 ] ] );
			if( dis_2 <= mesh.params.eps_2 )
			{
				mesh.tet_vertices[ v_id ].setFlag( eVertexFlags::Boundary );
				mesh.tet_vertices[ v_id ].setFlag( eVertexFlags::Cut, is_on_cut );
			}
		}
	};

	timer.start();
	std::vector<std::vector<std::pair<int, int>>> covered_fs_infos( input_faces.size() );
	for( int i = 0; i < track_surface_fs.size(); i++ )
	{
		if( mesh.tets[ i ].is_removed )
			continue;
		for( int j = 0; j < 4; j++ )
		{
			for( int f_id : track_surface_fs[ i ][ j ] )
				covered_fs_infos[ f_id ].push_back( std::make_pair( i, j ) );
		}
	}
	mesh.logger().logInfo( "time1 = %g", timer.getElapsedTime() );

	bool is_all_inserted = true;
	int cnt = 0;
	double time2 = 0;  // fortest
	double time3 = 0;
	double time4 = 0;
	double time5 = 0;
	double time6 = 0;
	for( int I = 0; I < b_edge_infos.size(); I++ )
	{
		const auto& e = b_edge_infos[ I ].first;
		auto& n_f_ids = b_edge_infos[ I ].second;  /// it is sorted
		bool is_on_cut = is_on_cut_edges[ I ];
		if( !is_again )
		{
			mesh.tet_vertices[ e[ 0 ] ].setFlag( eVertexFlags::Boundary );
			mesh.tet_vertices[ e[ 1 ] ].setFlag( eVertexFlags::Boundary );
			mesh.tet_vertices[ e[ 0 ] ].setFlag( eVertexFlags::Cut, is_on_cut );
			mesh.tet_vertices[ e[ 1 ] ].setFlag( eVertexFlags::Cut, is_on_cut );
		}
		//        b_edges.push_back(e);

		timer.start();
		/// double check neighbors
		for( int i = 0; i < n_f_ids.size(); i++ )
		{
			if( !is_face_inserted[ n_f_ids[ i ] ] )
			{
				n_f_ids.erase( n_f_ids.begin() + i );
				i--;
				break;
			}
		}
		time2 += timer.getElapsedTime();
		if( n_f_ids.empty() )
		{
			mesh.logger().logWarning( "FAIL n_f_ids.empty()" );
			continue;
		}

		timer.start();
		/// compute intersection
		std::vector<Vector3> points;
		std::map<std::array<int, 2>, int> map_edge_to_intersecting_point;
		std::vector<int> snapped_v_ids;
		std::vector<std::array<int, 3>> cut_fs;
		if( !insert_boundary_edges_get_intersecting_edges_and_points( covered_fs_infos, input_vertices, input_faces, e, n_f_ids, track_surface_fs, mesh, points,
			  map_edge_to_intersecting_point, snapped_v_ids, cut_fs, is_again ) )
		{
			for( int f_id : n_f_ids )
				is_face_inserted[ f_id ] = false;
			is_all_inserted = false;

			mesh.logger().logWarning( "FAIL insert_boundary_edges_get_intersecting_edges_and_points" );
			time3 += timer.getElapsedTime();
			continue;
		}
		time3 += timer.getElapsedTime();
		if( points.empty() )
		{  /// if all snapped
			record_boundary_info( points, snapped_v_ids, e, is_on_cut );
			cnt++;
			continue;
		}

		timer.start();
		/// subdivision
		std::vector<int> cut_t_ids;
		for( const auto& m : map_edge_to_intersecting_point )
		{
			const auto& tet_e = m.first;
			std::vector<int> tmp;
			setIntersection( mesh.tet_vertices[ tet_e[ 0 ] ].connTets, mesh.tet_vertices[ tet_e[ 1 ] ].connTets, tmp );
			cut_t_ids.insert( cut_t_ids.end(), tmp.begin(), tmp.end() );
		}
		vector_unique( cut_t_ids );
		std::vector<bool> is_mark_surface( cut_t_ids.size(), false );
		CutMesh empty_cut_mesh( mesh, Vector3( 0, 0, 0 ), std::array<Vector3, 3>() );
		//
		std::vector<MeshTet> new_tets;
		std::vector<std::array<std::vector<int>, 4>> new_track_surface_fs;
		std::vector<int> modified_t_ids;
		if( !subdivide_tets( -1, mesh, empty_cut_mesh, points, map_edge_to_intersecting_point, track_surface_fs, cut_t_ids, is_mark_surface, new_tets,
			  new_track_surface_fs, modified_t_ids ) )
		{
			bool is_inside_envelope = true;
			for( auto& f : cut_fs )
			{
				GEO2::index_t prev_facet = GEO2::NO_FACET;
				if( sample_triangle_and_check_is_out( { { mesh.tet_vertices[ f[ 0 ] ].pos, mesh.tet_vertices[ f[ 1 ] ].pos, mesh.tet_vertices[ f[ 2 ] ].pos } },
					  mesh.params.dd, mesh.params.eps_2, tree, prev_facet ) )
				{
					is_inside_envelope = false;
					break;
				}
			}
			if( !is_inside_envelope )
			{
				for( int f_id : n_f_ids )
					is_face_inserted[ f_id ] = false;
				for( auto& f : cut_fs )
				{
					mark_known_surface_fs( f, KNOWN_NOT_SURFACE );
				}
				known_not_surface_fs.insert( known_not_surface_fs.end(), cut_fs.begin(), cut_fs.end() );
				mesh.logger().logWarning( "FAIL subdivide_tets" );
			}
			else
			{
				for( auto& f : cut_fs )
					mark_known_surface_fs( f, KNOWN_SURFACE );
				known_surface_fs.insert( known_surface_fs.end(), cut_fs.begin(), cut_fs.end() );
				mesh.logger().logWarning( "SEMI-FAIL subdivide_tets" );
			}

			is_all_inserted = false;  // unless now
			time4 += timer.getElapsedTime();
			continue;
		}
		time4 += timer.getElapsedTime();

		//
		timer.start();
		push_new_tets( mesh, track_surface_fs, points, new_tets, new_track_surface_fs, modified_t_ids, is_again );
		time5 += timer.getElapsedTime();

		//
		/// mark boundary vertices
		/// notice, here we assume points list is inserted in the end of mesh.tet_vertices
		timer.start();
		record_boundary_info( points, snapped_v_ids, e, is_on_cut );
		time6 += timer.getElapsedTime();
		cnt++;
	}

	mesh.logger().logInfo( "uninsert boundary #e = %i/%zu", (int)b_edge_infos.size() - cnt, b_edge_infos.size() );
	mesh.logger().logInfo( "time2 = %g", time2 );
	mesh.logger().logInfo( "time3 = %g", time3 );
	mesh.logger().logInfo( "time4 = %g", time4 );
	mesh.logger().logInfo( "time5 = %g", time5 );
	mesh.logger().logInfo( "time6 = %g", time6 );
	return is_all_inserted;
}

bool floatTetWild::insert_boundary_edges_get_intersecting_edges_and_points( const std::vector<std::vector<std::pair<int, int>>>& covered_fs_infos,
  const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::array<int, 2>& e, const std::vector<int>& n_f_ids,
  std::vector<std::array<std::vector<int>, 4>>& track_surface_fs, Mesh& mesh, std::vector<Vector3>& points,
  std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point, std::vector<int>& snapped_v_ids, std::vector<std::array<int, 3>>& cut_fs, bool is_again )
{
	//    igl::Timer timer;

	auto is_cross = []( eOrientation a, eOrientation b )
	{
		if( ( a == eOrientation::Positive && b == eOrientation::Negative ) || ( a == eOrientation::Negative && b == eOrientation::Positive ) )
			return true;
		return false;
	};

	int t = get_t( input_vertices[ input_faces[ n_f_ids.front() ][ 0 ] ], input_vertices[ input_faces[ n_f_ids.front() ][ 1 ] ],
	  input_vertices[ input_faces[ n_f_ids.front() ][ 2 ] ] );
	std::array<Vector2, 2> evs_2d = { { to_2d( input_vertices[ e[ 0 ] ], t ), to_2d( input_vertices[ e[ 1 ] ], t ) } };
	Vector3 n = ( input_vertices[ input_faces[ n_f_ids.front() ][ 1 ] ] - input_vertices[ input_faces[ n_f_ids.front() ][ 0 ] ] )
				  .cross( input_vertices[ input_faces[ n_f_ids.front() ][ 2 ] ] - input_vertices[ input_faces[ n_f_ids.front() ][ 0 ] ] );
	n.normalize();
	const Vector3& pp = input_vertices[ input_faces[ n_f_ids.front() ][ 0 ] ];

	//    timer.start();
	std::vector<bool> is_visited( mesh.tets.size(), false );
	std::queue<int> t_ids_queue;
	/// find seed t_ids
	if( !is_again )
	{
		std::vector<int> t_ids;
		for( int f_id : n_f_ids )
		{
			for( const auto& info : covered_fs_infos[ f_id ] )
				t_ids.push_back( info.first );
		}
		for( int v_id : e )
			t_ids.insert( t_ids.end(), mesh.tet_vertices[ v_id ].connTets.begin(), mesh.tet_vertices[ v_id ].connTets.end() );
		vector_unique( t_ids );
		for( int t_id : t_ids )
		{
			t_ids_queue.push( t_id );
			is_visited[ t_id ] = true;
		}
	}
	else
	{
		__m256d min_e, max_e;
		get_bbox_face( input_vertices[ e[ 0 ] ], input_vertices[ e[ 0 ] ], input_vertices[ e[ 1 ] ], min_e, max_e );

		std::vector<int> t_ids;
		for( int f_id : n_f_ids )
		{
			for( const auto& info : covered_fs_infos[ f_id ] )
				t_ids.push_back( info.first );
		}
#ifdef FLOAT_TETWILD_USE_TBB
		tbb::concurrent_vector<int> tbb_t_ids;
		tbb::parallel_for( size_t( 0 ), mesh.tets.size(),
		  [ & ]( size_t t_id )
		  {
			  if( mesh.tets[ t_id ].is_removed )
				  return;
			  if( track_surface_fs[ t_id ][ 0 ].empty() && track_surface_fs[ t_id ][ 1 ].empty() && track_surface_fs[ t_id ][ 2 ].empty() &&
				  track_surface_fs[ t_id ][ 3 ].empty() )
				  return;

			  Vector3 min_t, max_t;
			  get_bbox_tet( mesh.tet_vertices[ mesh.tets[ t_id ][ 0 ] ].pos, mesh.tet_vertices[ mesh.tets[ t_id ][ 1 ] ].pos,
				mesh.tet_vertices[ mesh.tets[ t_id ][ 2 ] ].pos, mesh.tet_vertices[ mesh.tets[ t_id ][ 3 ] ].pos, min_t, max_t );
			  if( !is_bbox_intersected( min_e, max_e, min_t, max_t ) )
				  return;
			  tbb_t_ids.push_back( t_id );
		  } );
		t_ids.insert( t_ids.end(), tbb_t_ids.begin(), tbb_t_ids.end() );
#else
		for( int t_id = 0; t_id < mesh.tets.size(); t_id++ )
		{
			if( mesh.tets[ t_id ].is_removed )
				continue;
			if( track_surface_fs[ t_id ][ 0 ].empty() && track_surface_fs[ t_id ][ 1 ].empty() && track_surface_fs[ t_id ][ 2 ].empty() &&
				track_surface_fs[ t_id ][ 3 ].empty() )
				continue;

			__m256d min_t, max_t;
			get_bbox_tet( mesh.tet_vertices[ mesh.tets[ t_id ][ 0 ] ].pos, mesh.tet_vertices[ mesh.tets[ t_id ][ 1 ] ].pos,
			  mesh.tet_vertices[ mesh.tets[ t_id ][ 2 ] ].pos, mesh.tet_vertices[ mesh.tets[ t_id ][ 3 ] ].pos, min_t, max_t );
			if( !is_bbox_intersected( min_e, max_e, min_t, max_t ) )
				continue;

			t_ids.push_back( t_id );
		}
#endif

		vector_unique( t_ids );
		for( int t_id : t_ids )
		{
			t_ids_queue.push( t_id );
			is_visited[ t_id ] = true;
		}
	}
	//    time_e1+=timer.getElapsedTimeInSec();

	//    timer.start();
	std::vector<eOrientation> v_oris( mesh.tet_vertices.size(), eOrientation::Unknown );
	while( !t_ids_queue.empty() )
	{
		int t_id = t_ids_queue.front();
		t_ids_queue.pop();

		std::array<bool, 4> is_cut_vs = { { false, false, false, false } };
		for( int j = 0; j < 4; j++ )
		{
			/// check if contains
			if( track_surface_fs[ t_id ][ j ].empty() )
				continue;
			std::sort( track_surface_fs[ t_id ][ j ].begin(), track_surface_fs[ t_id ][ j ].end() );
			std::vector<int> tmp;
			std::set_intersection(
			  track_surface_fs[ t_id ][ j ].begin(), track_surface_fs[ t_id ][ j ].end(), n_f_ids.begin(), n_f_ids.end(), std::back_inserter( tmp ) );
			if( tmp.empty() )
				continue;

			/// check if cut through
			// check tri side of seg
			std::array<int, 3> f_v_ids = { { mesh.tets[ t_id ][ ( j + 1 ) % 4 ], mesh.tets[ t_id ][ ( j + 2 ) % 4 ], mesh.tets[ t_id ][ ( j + 3 ) % 4 ] } };
			std::array<Vector2, 3> fvs_2d = { { to_2d( mesh.tet_vertices[ f_v_ids[ 0 ] ].pos, n, pp, t ),
			  to_2d( mesh.tet_vertices[ f_v_ids[ 1 ] ].pos, n, pp, t ), to_2d( mesh.tet_vertices[ f_v_ids[ 2 ] ].pos, n, pp, t ) } };
			int cnt_pos = 0;
			int cnt_neg = 0;
			int cnt_on = 0;
			for( int k = 0; k < 3; k++ )
			{
				eOrientation& ori = v_oris[ f_v_ids[ k ] ];
				if( ori == eOrientation::Unknown )
					ori = Predicates::orient_2d( evs_2d[ 0 ], evs_2d[ 1 ], fvs_2d[ k ] );
				if( ori == eOrientation::Zero )
				{
					cnt_on++;
				}
				else
				{
					Scalar dis_2 = p_seg_squared_dist_3d( mesh.tet_vertices[ f_v_ids[ k ] ].pos, input_vertices[ e[ 0 ] ], input_vertices[ e[ 1 ] ] );
					if( dis_2 < mesh.params.eps_2_coplanar )
					{
						ori = eOrientation::Zero;
						cnt_on++;
						continue;
					}
					if( ori == eOrientation::Positive )
						cnt_pos++;
					else
						cnt_neg++;
				}
			}
			if( cnt_on >= 2 )
			{
				cut_fs.push_back( f_v_ids );
				std::sort( cut_fs.back().begin(), cut_fs.back().end() );
				continue;
			}
			if( cnt_neg == 0 || cnt_pos == 0 )
				continue;

			// check tri edge - seg intersection
			bool is_intersected = false;
			for( auto& p : evs_2d )
			{  /// first check if endpoints are contained inside the triangle
				if( is_p_inside_tri_2d( p, fvs_2d ) )
				{
					is_intersected = true;
					break;
				}
			}
			if( !is_intersected )
			{  /// then check if there's intersection
				for( int k = 0; k < 3; k++ )
				{
					// if cross
					if( !is_cross( v_oris[ f_v_ids[ k ] ], v_oris[ f_v_ids[ ( k + 1 ) % 3 ] ] ) )
						continue;
					// if already know intersect
					std::array<int, 2> tri_e = { { f_v_ids[ k ], f_v_ids[ ( k + 1 ) % 3 ] } };
					if( tri_e[ 0 ] > tri_e[ 1 ] )
						std::swap( tri_e[ 0 ], tri_e[ 1 ] );
					if( map_edge_to_intersecting_point.find( tri_e ) != map_edge_to_intersecting_point.end() )
					{
						is_intersected = true;
						break;
					}
					// if intersect
					double t2 = -1;
					if( seg_seg_intersection_2d( evs_2d, { { fvs_2d[ k ], fvs_2d[ ( k + 1 ) % 3 ] } }, t2 ) )
					{
						Vector3 p = ( 1 - t2 ) * mesh.tet_vertices[ f_v_ids[ k ] ].pos + t2 * mesh.tet_vertices[ f_v_ids[ ( k + 1 ) % 3 ] ].pos;
						double dis1 = ( p - mesh.tet_vertices[ f_v_ids[ k ] ].pos ).squaredNorm();
						double dis2 = ( p - mesh.tet_vertices[ f_v_ids[ ( k + 1 ) % 3 ] ].pos ).squaredNorm();
						//                        if (dis1 < SCALAR_ZERO_2) {
						if( dis1 < mesh.params.eps_2_coplanar )
						{
							v_oris[ f_v_ids[ k ] ] = eOrientation::Zero;
							is_intersected = true;
							break;
						}
						//                        if (dis2 < SCALAR_ZERO_2) {
						if( dis2 < mesh.params.eps_2_coplanar )
						{
							v_oris[ f_v_ids[ k ] ] = eOrientation::Zero;
							is_intersected = true;
							break;
						}
						points.push_back( p );
						map_edge_to_intersecting_point[ tri_e ] = points.size() - 1;
						is_intersected = true;
						break;
					}
				}
				if( !is_intersected )  /// no need to return false here
					continue;
			}

			std::sort( f_v_ids.begin(), f_v_ids.end() );
			cut_fs.push_back( f_v_ids );
			is_cut_vs[ ( j + 1 ) % 4 ] = true;
			is_cut_vs[ ( j + 2 ) % 4 ] = true;
			is_cut_vs[ ( j + 3 ) % 4 ] = true;
		}
		for( int j = 0; j < 4; j++ )
		{
			if( !is_cut_vs[ j ] )
				continue;
			for( int n_t_id : mesh.tet_vertices[ mesh.tets[ t_id ][ j ] ].connTets )
			{
				if( !is_visited[ n_t_id ] )
				{
					t_ids_queue.push( n_t_id );
					is_visited[ n_t_id ] = true;
				}
			}
		}
	}
	vector_unique( cut_fs );
	//    time_e2+=timer.getElapsedTimeInSec();

	//    timer.start();
	std::vector<std::array<int, 2>> tet_edges;
	for( const auto& f : cut_fs )
	{
		for( int j = 0; j < 3; j++ )
		{
			if( f[ j ] < f[ ( j + 1 ) % 3 ] )
				tet_edges.push_back( { { f[ j ], f[ ( j + 1 ) % 3 ] } } );
			else
				tet_edges.push_back( { { f[ ( j + 1 ) % 3 ], f[ j ] } } );
		}
	}
	vector_unique( tet_edges );
	//
	for( const auto& tet_e : tet_edges )
	{
		bool is_snapped = false;
		for( int j = 0; j < 2; j++ )
		{
			if( v_oris[ tet_e[ j ] ] == eOrientation::Zero )
			{
				snapped_v_ids.push_back( tet_e[ j ] );
				is_snapped = true;
			}
		}
		if( is_snapped )
			continue;
		if( map_edge_to_intersecting_point.find( tet_e ) != map_edge_to_intersecting_point.end() )
			continue;
		std::array<Vector2, 2> tri_evs_2d = {
		  { to_2d( mesh.tet_vertices[ tet_e[ 0 ] ].pos, n, pp, t ), to_2d( mesh.tet_vertices[ tet_e[ 1 ] ].pos, n, pp, t ) } };
		Scalar t_seg = -1;
		if( seg_line_intersection_2d( tri_evs_2d, evs_2d, t_seg ) )
		{
			points.push_back( ( 1 - t_seg ) * mesh.tet_vertices[ tet_e[ 0 ] ].pos + t_seg * mesh.tet_vertices[ tet_e[ 1 ] ].pos );
			map_edge_to_intersecting_point[ tet_e ] = points.size() - 1;
		}
	}
	vector_unique( snapped_v_ids );
	return true;
}

void floatTetWild::mark_surface_fs( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::vector<int>& input_tags,
  std::vector<std::array<std::vector<int>, 4>>& track_surface_fs, const std::vector<bool>& is_face_inserted,
  const std::vector<std::array<int, 3>>& known_surface_fs, const std::vector<std::array<int, 3>>& known_not_surface_fs,
  std::vector<std::array<int, 2>>& b_edges, Mesh& mesh, AABBWrapper& tree )
{
	auto is_on_bounded_side = [ & ]( const std::array<Vector2, 3>& ps_2d, const Vector2& c )
	{
		int cnt_pos = 0;
		int cnt_neg = 0;
		for( int j = 0; j < 3; j++ )
		{
			eOrientation ori = Predicates::orient_2d( ps_2d[ j ], ps_2d[ ( j + 1 ) % 3 ], c );
			if( ori == eOrientation::Positive )
				cnt_pos++;
			else if( ori == eOrientation::Negative )
				cnt_neg++;
		}

		if( cnt_neg > 0 && cnt_pos > 0 )
			return false;
		return true;
	};

	std::vector<std::array<bool, 4>> is_visited( track_surface_fs.size(), { { false, false, false, false } } );
	for( int t_id = 0; t_id < track_surface_fs.size(); t_id++ )
	{
		for( int j = 0; j < 4; j++ )
		{
			int ff_id = -1;
			int opp_t_id = -1;
			int k = -1;
			if( mesh.tets[ t_id ].is_surface_fs[ j ] == KNOWN_SURFACE || mesh.tets[ t_id ].is_surface_fs[ j ] == KNOWN_NOT_SURFACE )
			{
				opp_t_id = get_opp_t_id( t_id, j, mesh );
				if( opp_t_id < 0 )
				{
					mesh.tets[ t_id ].is_surface_fs[ j ] = NOT_SURFACE;
					continue;
				}
				k =
				  get_local_f_id( opp_t_id, mesh.tets[ t_id ][ ( j + 1 ) % 4 ], mesh.tets[ t_id ][ ( j + 2 ) % 4 ], mesh.tets[ t_id ][ ( j + 3 ) % 4 ], mesh );
				is_visited[ t_id ][ j ] = true;
				is_visited[ opp_t_id ][ k ] = true;
				if( mesh.tets[ t_id ].is_surface_fs[ j ] == KNOWN_NOT_SURFACE || track_surface_fs[ t_id ][ j ].empty() )
				{
					mesh.tets[ t_id ].is_surface_fs[ j ] = NOT_SURFACE;
					mesh.tets[ opp_t_id ].is_surface_fs[ k ] = NOT_SURFACE;
					continue;
				}
				else
					ff_id = track_surface_fs[ t_id ][ j ].front();
			}
			else
			{
				if( mesh.tets[ t_id ].is_surface_fs[ j ] != NOT_SURFACE || is_visited[ t_id ][ j ] )
					continue;
				is_visited[ t_id ][ j ] = true;
				if( track_surface_fs[ t_id ][ j ].empty() )
					continue;

				auto& f_ids = track_surface_fs[ t_id ][ j ];

				auto& tp1_3d = mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 1 ) % 4 ] ].pos;
				auto& tp2_3d = mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 2 ) % 4 ] ].pos;
				auto& tp3_3d = mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 3 ) % 4 ] ].pos;

				double eps_2 = ( mesh.params.eps + mesh.params.eps_simplification ) / 2;
				double dd = ( mesh.params.dd + mesh.params.dd_simplification ) / 2;
				eps_2 *= eps_2;
				GEO2::index_t prev_facet = GEO2::NO_FACET;
				if( sample_triangle_and_check_is_out( { { tp1_3d, tp2_3d, tp3_3d } }, dd, eps_2, tree, prev_facet ) )
					continue;
				else
					ff_id = track_surface_fs[ t_id ][ j ].front();

				opp_t_id = get_opp_t_id( t_id, j, mesh );
				if( opp_t_id < 0 )
				{
					mesh.tets[ t_id ].is_surface_fs[ j ] = NOT_SURFACE;
					continue;
				}
				k =
				  get_local_f_id( opp_t_id, mesh.tets[ t_id ][ ( j + 1 ) % 4 ], mesh.tets[ t_id ][ ( j + 2 ) % 4 ], mesh.tets[ t_id ][ ( j + 3 ) % 4 ], mesh );
				is_visited[ opp_t_id ][ k ] = true;
			}

			myassert( ff_id >= 0, "ff_id<0!!!" );  // fortest

			mesh.tets[ t_id ].surface_tags[ j ] = input_tags[ ff_id ];
			mesh.tets[ opp_t_id ].surface_tags[ k ] = input_tags[ ff_id ];

			auto& fv1 = input_vertices[ input_faces[ ff_id ][ 0 ] ];
			auto& fv2 = input_vertices[ input_faces[ ff_id ][ 1 ] ];
			auto& fv3 = input_vertices[ input_faces[ ff_id ][ 2 ] ];
			//
			eOrientation ori = Predicates::orient_3d( fv1, fv2, fv3, mesh.tet_vertices[ mesh.tets[ t_id ][ j ] ].pos );
			eOrientation opp_ori = Predicates::orient_3d( fv1, fv2, fv3, mesh.tet_vertices[ mesh.tets[ opp_t_id ][ k ] ].pos );
			//
			if( ori == eOrientation::Positive && opp_ori == eOrientation::Negative || ori == eOrientation::Negative && opp_ori == eOrientation::Positive )
			{
				mesh.tets[ t_id ].is_surface_fs[ j ] = (int8_t)ori;
				mesh.tets[ opp_t_id ].is_surface_fs[ k ] = (int8_t)opp_ori;
				continue;
			}
			//
			if( ori == eOrientation::Zero && opp_ori != eOrientation::Zero )
			{
				mesh.tets[ t_id ].is_surface_fs[ j ] = -(int8_t)opp_ori;
				mesh.tets[ opp_t_id ].is_surface_fs[ k ] = (int8_t)opp_ori;
				continue;
			}
			if( opp_ori == eOrientation::Zero && ori != eOrientation::Zero )
			{
				mesh.tets[ t_id ].is_surface_fs[ j ] = (int8_t)ori;
				mesh.tets[ opp_t_id ].is_surface_fs[ k ] = -(int8_t)ori;
				continue;
			}
			//
			if( ori == opp_ori )
			{
				Vector3 n = ( fv2 - fv1 ).cross( fv3 - fv1 );
				n.normalize();
				Scalar dist = n.dot( mesh.tet_vertices[ mesh.tets[ t_id ][ j ] ].pos - fv1 );
				Scalar opp_dist = n.dot( mesh.tet_vertices[ mesh.tets[ opp_t_id ][ k ] ].pos - fv1 );
				if( ori == eOrientation::Zero )
				{
					auto& tv1 = mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 1 ) % 4 ] ].pos;
					auto& tv2 = mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 2 ) % 4 ] ].pos;
					auto& tv3 = mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 3 ) % 4 ] ].pos;
					Vector3 nt;
					if( Predicates::orient_3d( tv1, tv2, tv3, mesh.tet_vertices[ mesh.tets[ t_id ][ j ] ].pos ) == eOrientation::Positive )
						nt = ( tv2 - tv1 ).cross( tv3 - tv1 );
					else
						nt = ( tv3 - tv1 ).cross( tv2 - tv1 );
					nt.normalize();
					if( n.dot( nt ) > 0 )
					{
						mesh.tets[ opp_t_id ].is_surface_fs[ k ] = 1;
						mesh.tets[ t_id ].is_surface_fs[ j ] = -1;
					}
					else
					{
						mesh.tets[ opp_t_id ].is_surface_fs[ k ] = -1;
						mesh.tets[ t_id ].is_surface_fs[ j ] = 1;
					}
				}
				else
				{
					if( dist < opp_dist )
					{
						mesh.tets[ opp_t_id ].is_surface_fs[ k ] = (int8_t)opp_ori;
						mesh.tets[ t_id ].is_surface_fs[ j ] = -(int8_t)ori;
					}
					else
					{
						mesh.tets[ opp_t_id ].is_surface_fs[ k ] = -(int8_t)opp_ori;
						mesh.tets[ t_id ].is_surface_fs[ j ] = (int8_t)ori;
					}
				}
			}
			//
		}
	}

	mesh.logger().logDebug( "known_surface_fs.size = %zu, known_not_surface_fs.size = %zu", known_surface_fs.size(), known_not_surface_fs.size() );
	if( known_surface_fs.empty() && known_not_surface_fs.empty() )
		return;

	// b_edges
	std::vector<std::array<int, 2>> tmp_edges;
	for( const auto& f : known_surface_fs )
	{
		for( int j = 0; j < 3; j++ )
		{
			if( f[ j ] < f[ ( j + 1 ) % 3 ] )
				tmp_edges.push_back( { { f[ j ], f[ ( j + 1 ) % 3 ] } } );
			else
				tmp_edges.push_back( { { f[ ( j + 1 ) % 3 ], f[ j ] } } );
		}
	}
	for( const auto& f : known_not_surface_fs )
	{
		for( int j = 0; j < 3; j++ )
		{
			if( f[ j ] < f[ ( j + 1 ) % 3 ] )
				tmp_edges.push_back( { { f[ j ], f[ ( j + 1 ) % 3 ] } } );
			else
				tmp_edges.push_back( { { f[ ( j + 1 ) % 3 ], f[ j ] } } );
		}
	}
	vector_unique( tmp_edges );

	for( auto& e : tmp_edges )
	{
		int cnt = 0;
		std::vector<int> n_t_ids;
		setIntersection( mesh.tet_vertices[ e[ 0 ] ].connTets, mesh.tet_vertices[ e[ 1 ] ].connTets, n_t_ids );
		for( int t_id : n_t_ids )
		{
			for( int j = 0; j < 4; j++ )
			{
				if( mesh.tets[ t_id ][ j ] != e[ 0 ] && mesh.tets[ t_id ][ j ] != e[ 1 ] )
				{
					if( mesh.tets[ t_id ].is_surface_fs[ j ] != NOT_SURFACE )
					{
						cnt++;
						if( cnt > 2 )
							break;
					}
				}
				if( cnt > 2 )
					break;
			}
		}
		//        cout<<"cnt = "<<cnt<<endl;
		if( cnt == 2 )
		{
			b_edges.push_back( e );
			mesh.tet_vertices[ e[ 0 ] ].setFlag( eVertexFlags::Boundary );
			//            cout<<"b_edges.push_back(e);"<<endl;
		}
	}
}

bool floatTetWild::is_uninserted_face_covered(
  int uninserted_f_id, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::vector<int>& cut_t_ids, Mesh& mesh )
{
	std::array<Vector3, 3> vs = { { input_vertices[ input_faces[ uninserted_f_id ][ 0 ] ], input_vertices[ input_faces[ uninserted_f_id ][ 1 ] ],
	  input_vertices[ input_faces[ uninserted_f_id ][ 2 ] ] } };
	std::vector<GEO2::vec3> ps;
	sample_triangle( vs, ps, mesh.params.dd );

	std::vector<int> n_t_ids;
	for( int t_id : cut_t_ids )
	{
		for( int j = 0; j < 4; j++ )
			n_t_ids.insert(
			  n_t_ids.end(), mesh.tet_vertices[ mesh.tets[ t_id ][ j ] ].connTets.begin(), mesh.tet_vertices[ mesh.tets[ t_id ][ j ] ].connTets.end() );
	}
	vector_unique( n_t_ids );

	std::vector<std::array<int, 3>> faces;
	for( int t_id : n_t_ids )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( mesh.tets[ t_id ].is_surface_fs[ j ] != NOT_SURFACE )
			{
				faces.push_back( { { mesh.tets[ t_id ][ ( j + 1 ) % 4 ], mesh.tets[ t_id ][ ( j + 2 ) % 4 ], mesh.tets[ t_id ][ ( j + 3 ) % 4 ] } } );
				std::sort( faces.back().begin(), faces.back().end() );
			}
		}
	}
	vector_unique( faces );

	for( auto& p : ps )
	{
		bool is_valid = false;
		for( auto& f : faces )
		{
			double dis_2 = GEO2::point_triangle_squared_distance(
			  p, to_geo_p( mesh.tet_vertices[ f[ 0 ] ].pos ), to_geo_p( mesh.tet_vertices[ f[ 1 ] ].pos ), to_geo_p( mesh.tet_vertices[ f[ 2 ] ].pos ) );
			if( dis_2 < mesh.params.eps_2 )
			{
				is_valid = true;
				break;
			}
		}
		if( !is_valid )
			return false;
	}

	mesh.logger().logDebug( "covered!!!!!!!" );
	return true;
}

int floatTetWild::get_opp_t_id( int t_id, int j, const Mesh& mesh )
{
	std::vector<int> tmp;
	setIntersection( mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 1 ) % 4 ] ].connTets, mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 2 ) % 4 ] ].connTets,
	  mesh.tet_vertices[ mesh.tets[ t_id ][ ( j + 3 ) % 4 ] ].connTets, tmp );
	if( tmp.size() == 2 )
		return tmp[ 0 ] == t_id ? tmp[ 1 ] : tmp[ 0 ];
	else
		return -1;
}

void floatTetWild::myassert( bool b, const std::string& s )
{
	if( b == false )
	{
		throw std::logic_error( s );
	}
}

floatTetWild::eOrientation floatTetWild::orient_rational( const Vector3_r& p1, const Vector3_r& p2, const Vector3_r& p3, const Vector3_r& p )
{
	auto nv = ( p2 - p1 ).cross( p3 - p1 );
	triwild::Rational res = nv.dot( p - p1 );
	if( res == 0 )
		return eOrientation::Zero;
	if( res < 0 )
		return eOrientation::Positive;
	else
		return eOrientation::Negative;
}