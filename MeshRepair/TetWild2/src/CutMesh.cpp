// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

//
// Created by Yixin Hu on 9/12/19.
//

#include "stdafx.h"
#include "CutMesh.h"
#include "TriangleInsertion.h"
#include "LocalOperations.h"
#include "../external/Predicates.h"
#include "intersections.h"
#include <igl/Timer.h>

void floatTetWild::CutMesh::construct( const std::vector<int>& cut_t_ids )
{
	v_ids.clear();
	for( int t_id : cut_t_ids )
		for( int j = 0; j < 4; j++ )
			v_ids.addSorted( mesh.tets[ t_id ][ j ] );

	for( int i = 0; i < v_ids.size(); i++ )
		map_v_ids.setAt( v_ids[ i ], i );

	tets.resize( cut_t_ids.size() );
	for( int i = 0; i < cut_t_ids.size(); i++ )
	{
		tets[ i ] = { { map_v_ids[ mesh.tets[ cut_t_ids[ i ] ][ 0 ] ], map_v_ids[ mesh.tets[ cut_t_ids[ i ] ][ 1 ] ],
		  map_v_ids[ mesh.tets[ cut_t_ids[ i ] ][ 2 ] ], map_v_ids[ mesh.tets[ cut_t_ids[ i ] ][ 3 ] ] } };
	}
}

bool floatTetWild::CutMesh::snap_to_plane()
{
	bool snapped = false;
	to_plane_dists.resize( map_v_ids.size() );
	is_snapped.resize( map_v_ids.size() );
	for( auto& v : map_v_ids )
	{
		int v_id = v.first;
		int lv_id = v.second;

		const eOrientation ori = Predicates::orient_3d( p_vs[ 0 ], p_vs[ 1 ], p_vs[ 2 ], mesh.tet_vertices[ v_id ].pos );
		if( ori == eOrientation::Zero )
		{
			to_plane_dists[ lv_id ] = 0;
			continue;
		}
		to_plane_dists[ lv_id ] = get_to_plane_dist( mesh.tet_vertices[ v_id ].pos );
		if( ori == eOrientation::Positive && to_plane_dists[ lv_id ] > 0 || ori == eOrientation::Negative && to_plane_dists[ lv_id ] < 0 )
		{
			//            cout<<"reverted!!! "<<to_plane_dists[lv_id]<<endl;
			to_plane_dists[ lv_id ] = -to_plane_dists[ lv_id ];
		}

		if( std::fabs( to_plane_dists[ lv_id ] ) < mesh.params.eps_coplanar )
		{
			is_snapped[ lv_id ] = true;
			snapped = true;
		}
	}

	revert_totally_snapped_tets( 0, tets.size() );

	return snapped;
}

void floatTetWild::CutMesh::expand( std::vector<int>& cut_t_ids )
{
	std::vector<std::vector<int>> conn_tets( v_ids.size() );
	for( int i = 0; i < tets.size(); i++ )
	{
		for( int j = 0; j < 4; j++ )
			conn_tets[ tets[ i ][ j ] ].push_back( i );
	}

	std::vector<std::array<int, 4>> opp_t_ids( tets.size(), { { -1, -1, -1, -1 } } );
	for( int i = 0; i < tets.size(); i++ )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( opp_t_ids[ i ][ j ] >= 0 )
				continue;

			std::vector<int> n_t_ids;
			set_intersection_sorted(
			  conn_tets[ tets[ i ][ ( j + 1 ) % 4 ] ], conn_tets[ tets[ i ][ ( j + 2 ) % 4 ] ], conn_tets[ tets[ i ][ ( j + 3 ) % 4 ] ], n_t_ids );

			assert( !n_t_ids.empty() );
			if( n_t_ids.size() < 2 )
				continue;

			int n_t_id = n_t_ids[ 0 ] == i ? n_t_ids[ 1 ] : n_t_ids[ 0 ];
			opp_t_ids[ i ][ j ] = n_t_id;
			for( int k = 0; k < 4; k++ )
			{
				if( tets[ n_t_id ][ k ] != tets[ i ][ ( j + 1 ) % 4 ] && tets[ n_t_id ][ k ] != tets[ i ][ ( j + 2 ) % 4 ] &&
					tets[ n_t_id ][ k ] != tets[ i ][ ( j + 3 ) % 4 ] )
				{
					opp_t_ids[ n_t_id ][ k ] = i;
					break;
				}
			}
		}
	}

	while( true )
	{
		std::vector<std::array<int, 5>> new_opp_t_ids;
		for( int i = 0; i < opp_t_ids.size(); i++ )
		{
			for( int j = 0; j < 4; j++ )
			{
				if( opp_t_ids[ i ][ j ] >= 0 )
					continue;
				if( !is_snapped[ tets[ i ][ ( j + 1 ) % 4 ] ] && !is_snapped[ tets[ i ][ ( j + 2 ) % 4 ] ] && !is_snapped[ tets[ i ][ ( j + 3 ) % 4 ] ] )
					continue;

				int n_gt_id = get_opp_t_id( cut_t_ids[ i ], j, mesh );
				if( n_gt_id < 0 )
					continue;
				new_opp_t_ids.push_back( { { -1, -1, -1, -1, n_gt_id } } );
				for( int k = 0; k < 4; k++ )
				{
					if( mesh.tets[ n_gt_id ][ k ] != v_ids[ tets[ i ][ ( j + 1 ) % 4 ] ] && mesh.tets[ n_gt_id ][ k ] != v_ids[ tets[ i ][ ( j + 2 ) % 4 ] ] &&
						mesh.tets[ n_gt_id ][ k ] != v_ids[ tets[ i ][ ( j + 3 ) % 4 ] ] )
					{
						new_opp_t_ids.back()[ k ] = i;
						break;
					}
				}
			}
		}
		if( new_opp_t_ids.empty() )
			return;

		std::sort(
		  new_opp_t_ids.begin(), new_opp_t_ids.end(), [ & ]( const std::array<int, 5>& a, const std::array<int, 5>& b ) { return a.back() < b.back(); } );
		for( int i = 0; i < new_opp_t_ids.size() - 1; i++ )
		{
			if( new_opp_t_ids[ i ].back() == new_opp_t_ids[ i + 1 ].back() )
			{
				for( int j = 0; j < 4; j++ )
				{
					if( new_opp_t_ids[ i ][ j ] >= 0 )
					{
						new_opp_t_ids[ i + 1 ][ j ] = new_opp_t_ids[ i ][ j ];
						break;
					}
				}
				new_opp_t_ids.erase( new_opp_t_ids.begin() + i );
				i--;
			}
		}

		const int old_tets_size = tets.size();
		for( int i = 0; i < new_opp_t_ids.size(); i++ )
		{
			///
			int cnt_on = 0;
			for( int j = 0; j < 4; j++ )
			{
				int v_id = mesh.tets[ new_opp_t_ids[ i ].back() ][ j ];
				int foundVert;
				if( map_v_ids.tryLookup( v_id, foundVert ) && is_v_on_plane( foundVert ) )
					cnt_on++;
			}
			if( cnt_on != 3 )
			{
				int cnt_pos = 0;
				int cnt_neg = 0;
				for( int j = 0; j < 4; j++ )
				{
					const eOrientation ori =
					  Predicates::orient_3d( p_vs[ 0 ], p_vs[ 1 ], p_vs[ 2 ], mesh.tet_vertices[ mesh.tets[ new_opp_t_ids[ i ].back() ][ j ] ].pos );
					if( ori == eOrientation::Positive )
						cnt_pos++;
					else if( ori == eOrientation::Negative )
						cnt_neg++;
				}
				if( cnt_neg == 0 || cnt_pos == 0 )
					continue;
			}

			///
			cut_t_ids.push_back( new_opp_t_ids[ i ].back() );

			int t_id = tets.size();
			tets.emplace_back();
			auto& t = tets.back();
			const int old_v_ids_size = v_ids.size();
			for( int j = 0; j < 4; j++ )
			{
				const int v_id = mesh.tets[ new_opp_t_ids[ i ].back() ][ j ];
				int lv_id;
				auto placed = map_v_ids.emplace( v_id );
				if( placed.second )
				{
					lv_id = (int)v_ids.size();
					v_ids.addUnsorted( v_id );
					*placed.first = lv_id;
					double dist = get_to_plane_dist( mesh.tet_vertices[ v_id ].pos );
					const eOrientation ori = Predicates::orient_3d( p_vs[ 0 ], p_vs[ 1 ], p_vs[ 2 ], mesh.tet_vertices[ v_id ].pos );
					if( ( ori == eOrientation::Negative && dist < 0 )  // todo: change get_to_plane_dist return value sign
						|| ( ori == eOrientation::Positive && dist > 0 ) )
						dist = -dist;
					else if( ori == eOrientation::Zero )
						dist = 0;
					to_plane_dists.push_back( dist );

					if( ori != eOrientation::Zero && std::fabs( to_plane_dists[ lv_id ] ) < mesh.params.eps_coplanar )
						is_snapped.push_back( true );
					else
						is_snapped.push_back( false );
					conn_tets.emplace_back();
				}
				else
					lv_id = *placed.first;
				t[ j ] = lv_id;
				conn_tets[ lv_id ].push_back( t_id );
			}

			opp_t_ids.push_back( { { new_opp_t_ids[ i ][ 0 ], new_opp_t_ids[ i ][ 1 ], new_opp_t_ids[ i ][ 2 ], new_opp_t_ids[ i ][ 3 ] } } );
			for( int j = 0; j < 4; j++ )
			{
				if( opp_t_ids.back()[ j ] < 0 )
				{
					if( t[ ( j + 1 ) % 4 ] < old_v_ids_size && t[ ( j + 2 ) % 4 ] < old_v_ids_size && t[ ( j + 3 ) % 4 ] < old_v_ids_size )
					{
						std::vector<int> tmp;
						set_intersection( conn_tets[ t[ ( j + 1 ) % 4 ] ], conn_tets[ t[ ( j + 2 ) % 4 ] ], conn_tets[ t[ ( j + 3 ) % 4 ] ], tmp );
						if( tmp.size() == 1 )  //一个三角形刚好填了一个v字的缺口，所有点都在v_ids中，但三角形不在cut_mesh中
							continue;
						opp_t_ids.back()[ j ] = tmp[ 0 ] == t_id ? tmp[ 1 ] : tmp[ 0 ];
					}
					else
						continue;
				}
				int opp_t_id = opp_t_ids.back()[ j ];
				for( int k = 0; k < 4; k++ )
				{
					if( tets[ opp_t_id ][ k ] != t[ ( j + 1 ) % 4 ] && tets[ opp_t_id ][ k ] != t[ ( j + 2 ) % 4 ] &&
						tets[ opp_t_id ][ k ] != t[ ( j + 3 ) % 4 ] )
					{
						opp_t_ids[ opp_t_id ][ k ] = t_id;
						break;
					}
				}
			}
		}
		if( old_tets_size == tets.size() )
			break;

		revert_totally_snapped_tets( old_tets_size, tets.size() );
	}
}

void floatTetWild::CutMesh::expand_new( std::vector<int>& cut_t_ids, size_t countTets )
{
	const int t = get_t( p_vs[ 0 ], p_vs[ 1 ], p_vs[ 2 ] );
	const std::array<Vector2, 3> tri_2d = { { to_2d( p_vs[ 0 ], t ), to_2d( p_vs[ 1 ], t ), to_2d( p_vs[ 2 ], t ) } };

	std::vector<bool> is_in_cutmesh( countTets, false );
	for( int t_id : cut_t_ids )
		is_in_cutmesh[ t_id ] = true;

	int cnt_loop = 0;
	std::vector<bool> is_interior( v_ids.size(), false );
	while( true )
	{
		cnt_loop++;
		/////
		std::vector<bool> is_visited( countTets, false );
		for( int t_id : cut_t_ids )
			is_visited[ t_id ] = true;

		/////
		int old_cut_t_ids = cut_t_ids.size();
		for( auto m : map_v_ids )
		{
			int gv_id = m.first;
			int lv_id = m.second;

			if( is_interior[ lv_id ] )
				continue;
			if( !is_snapped[ lv_id ] )
				continue;

			bool is_in = true;
			for( int gt_id : mesh.tet_vertices[ gv_id ].connTets )
			{
				if( is_in_cutmesh[ gt_id ] )
					continue;
				is_in = false;

				if( is_visited[ gt_id ] )
					continue;
				is_visited[ gt_id ] = true;

				///
				int cnt = 0;
				int cnt_on = 0;
				for( int j = 0; j < 4; j++ )
				{
					int tmp_gv_id = mesh.tets[ gt_id ][ j ];
					int found;
					if( map_v_ids.tryLookup( tmp_gv_id, found ) )
					{
						cnt++;
						if( is_v_on_plane( found ) )
							cnt_on++;
					}
				}
				if( cnt < 3 )
					continue;
				int cnt_pos = 0;
				int cnt_neg = 0;
				for( int j = 0; j < 4; j++ )
				{
					eOrientation ori = Predicates::orient_3d( p_vs[ 0 ], p_vs[ 1 ], p_vs[ 2 ], mesh.tet_vertices[ mesh.tets[ gt_id ][ j ] ].pos );
					if( ori == eOrientation::Positive )
						cnt_pos++;
					else if( ori == eOrientation::Negative )
						cnt_neg++;
				}
				if( cnt_neg == 0 || cnt_pos == 0 )
					continue;

				bool is_overlapped = false;
				std::array<Vector2, 4> tet_2d;
				for( int j = 0; j < 4; j++ )
				{
					Scalar dist = get_to_plane_dist( mesh.tet_vertices[ mesh.tets[ gt_id ][ j ] ].pos );
					Vector3 proj_p = mesh.tet_vertices[ mesh.tets[ gt_id ][ j ] ].pos - dist * p_n;
					tet_2d[ j ] = to_2d( proj_p, t );
				}
				for( int j = 0; j < 4; j++ )
				{
					if( is_tri_tri_cutted_2d( { { tet_2d[ ( j + 1 ) % 4 ], tet_2d[ ( j + 2 ) % 4 ], tet_2d[ ( j + 3 ) % 4 ] } }, tri_2d ) )
					{
						is_overlapped = true;
						break;
					}
				}
				if( !is_overlapped )
					continue;

				///
				cut_t_ids.push_back( gt_id );
				is_in_cutmesh[ gt_id ] = true;

				///
				tets.emplace_back();
				auto& t = tets.back();
				for( int j = 0; j < 4; j++ )
				{
					int new_gv_id = mesh.tets[ gt_id ][ j ];
					int new_lv_id;
					auto placement = map_v_ids.emplace( new_gv_id );
					if( placement.second )
					{
						new_lv_id = v_ids.size();
						*placement.first = new_lv_id;

						v_ids.addUnsorted( new_gv_id );
						is_interior.push_back( false );
						//
						double dist = get_to_plane_dist( mesh.tet_vertices[ new_gv_id ].pos );
						const eOrientation ori = Predicates::orient_3d( p_vs[ 0 ], p_vs[ 1 ], p_vs[ 2 ], mesh.tet_vertices[ new_gv_id ].pos );
						if( ( ori == eOrientation::Negative && dist < 0 ) || ( ori == eOrientation::Positive && dist > 0 ) )
							dist = -dist;
						else if( ori == eOrientation::Zero )
							dist = 0;
						to_plane_dists.push_back( dist );
						//
						if( ori != eOrientation::Zero && std::fabs( to_plane_dists[ new_lv_id ] ) < mesh.params.eps_coplanar )
							is_snapped.push_back( true );
						else
							is_snapped.push_back( false );
						is_projected.push_back( false );
					}
					else
						new_lv_id = *placement.first;
					t[ j ] = new_lv_id;
				}
			}
			if( is_in )
				is_interior[ lv_id ] = true;
		}
		if( cut_t_ids.size() == old_cut_t_ids )
			break;
	}
	revert_totally_snapped_tets( 0, tets.size() );
}

int floatTetWild::CutMesh::project_to_plane( int input_vertices_size )
{
	is_projected.resize( v_ids.size(), false );

	int cnt = 0;
	for( int i = 0; i < is_snapped.size(); i++ )
	{
		if( !is_snapped[ i ] || is_projected[ i ] )
			continue;
		if( v_ids[ i ] < input_vertices_size )
			continue;
		Scalar dist = get_to_plane_dist( mesh.tet_vertices[ v_ids[ i ] ].pos );
		Vector3 proj_p = mesh.tet_vertices[ v_ids[ i ] ].pos - p_n * dist;
		bool is_snappable = true;
		for( int t_id : mesh.tet_vertices[ v_ids[ i ] ].connTets )
		{
			int j = mesh.tets[ t_id ].find( v_ids[ i ] );
			if( is_inverted( mesh, t_id, j, proj_p ) )
			{
				is_snappable = false;
				break;
			}
		}
		if( is_snappable )
		{
			mesh.tet_vertices[ v_ids[ i ] ].pos = proj_p;
			is_projected[ i ] = true;
			to_plane_dists[ i ] = get_to_plane_dist( proj_p );
			cnt++;
		}
	}
	return cnt;
}

void floatTetWild::CutMesh::revert_totally_snapped_tets( int a, int b )
{
	int cnt = 0;
	for( int i = a; i < b; i++ )
	{
		const auto& t = tets[ i ];
		if( is_v_on_plane( t[ 0 ] ) && is_v_on_plane( t[ 1 ] ) && is_v_on_plane( t[ 2 ] ) && is_v_on_plane( t[ 3 ] ) )
		{
			auto tmp_t = t;
			std::sort( tmp_t.begin(), tmp_t.end(), [ & ]( int a, int b ) { return fabs( to_plane_dists[ a ] ) > fabs( to_plane_dists[ b ] ); } );
			for( int j = 0; j < 3; j++ )
			{
				if( is_snapped[ tmp_t[ j ] ] == true )
				{
					is_snapped[ tmp_t[ j ] ] = false;
					break;
				}
			}
		}
	}
}

bool floatTetWild::CutMesh::get_intersecting_edges_and_points(
  std::vector<Vector3>& points, FlatEdgeMap& map_edge_to_intersecting_point, std::vector<int>& subdivide_t_ids )
{
	tmpEdges.clear();
	for( auto& t : tets )
	{
		for( int j = 0; j < 3; j++ )
		{
			tmpEdges.add( t[ 0 ], t[ j + 1 ] );
			tmpEdges.add( t[ j + 1 ], t[ mod3( j + 1 ) + 1 ] );
		}
	}

	std::vector<int> e_v_ids;
	for( int i = 0; i < tmpEdges.size(); i++ )
	{
		const auto& e = tmpEdges[ i ];
		if( is_v_on_plane( e[ 0 ] ) || is_v_on_plane( e[ 1 ] ) )
			continue;
		if( to_plane_dists[ e[ 0 ] ] * to_plane_dists[ e[ 1 ] ] >= 0 )
			continue;

		int v1_id = v_ids[ e[ 0 ] ];
		int v2_id = v_ids[ e[ 1 ] ];
		Vector3 p;
		Scalar _;
		bool is_result = seg_plane_intersection( mesh.tet_vertices[ v1_id ].pos, mesh.tet_vertices[ v2_id ].pos, p_vs[ 0 ], p_n, p, _ );
		if( !is_result )
		{
			mesh.logger().logWarning( "seg_plane_intersection no result!" );
			return false;
		}

		map_edge_to_intersecting_point.setAt( v1_id, v2_id, (int)points.size() );
		points.push_back( p );
		e_v_ids.push_back( v1_id );
		e_v_ids.push_back( v2_id );
	}
	vector_unique( e_v_ids );
	for( int v_id : e_v_ids )
		subdivide_t_ids.insert( subdivide_t_ids.end(), mesh.tet_vertices[ v_id ].connTets.begin(), mesh.tet_vertices[ v_id ].connTets.end() );
	vector_unique( subdivide_t_ids );

	return true;
}

bool floatTetWild::CutMesh::check()
{
	return true;
	/*
	bool is_good = true;
	for( auto& m : map_v_ids )
	{
		int gv_id = m.first;
		int lv_id = m.second;

		Scalar dist = get_to_plane_dist( mesh.tet_vertices[ gv_id ].pos );
		if( std::fabs( dist ) < mesh.params.eps_coplanar && to_plane_dists[ lv_id ] != 0 && !is_snapped[ lv_id ] )
		{
			cout << "wrong vertex in cut mesh" << endl;
			cout << dist << endl;
			cout << to_plane_dists[ lv_id ] << endl;
			cout << is_snapped[ lv_id ] << endl;
			is_good = false;
		}
	}
	return is_good; */
}