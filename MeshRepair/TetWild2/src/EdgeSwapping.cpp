// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
#include "stdafx.h"
#include "EdgeSwapping.h"
#include "LocalOperations.h"
#include "MeshImprovement.h"
#include "SmallBuffer.h"
#include "Bitset32.h"

namespace floatTetWild
{
	bool is_es_check = false;
}

void floatTetWild::edge_swapping( Mesh& mesh )
{
	auto& tet_vertices = mesh.tet_vertices;
	auto& tets = mesh.tets;

	int counter = 0;
	int suc_counter3 = 0;
	int suc_counter4 = 0;
	int suc_counter5 = 0;

	mesh.reset_t_empty_start();
	mesh.reset_v_empty_start();

	auto is_swappable = [ & ]( int v1_id, int v2_id, const std::vector<int>& n12_t_ids )
	{
		if( n12_t_ids.size() < 3 || n12_t_ids.size() > 5 )
			return false;
		if( !is_valid_edge( mesh, v1_id, v2_id, n12_t_ids ) )
			return false;
		if( is_surface_edge( mesh, v1_id, v2_id, n12_t_ids ) )
			return false;
		if( is_bbox_edge( mesh, v1_id, v2_id, n12_t_ids ) )
			return false;
		return true;
	};

	EdgeSwappingBuffers& buffers = mesh.edgeSwappingBuffers;
	// Apparently, edge_swapping doesn't run in parallel with findNewPosBuffers; reusing an old buffer of the same type
	EdgesSet& edges = mesh.findNewPosBuffers[ 0 ].edgesTemp;
	get_all_edges( mesh, edges );

	std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_l>& es_queue = buffers.es_queue;
	edges.enumerate( [ & ]( int e0, int e1 ) { es_queue.push( ElementInQueue( e0, e1, get_edge_length_2( mesh, e0, e1 ) ) ); } );
	edges.clear();

	std::vector<int>& n12_t_ids = buffers.n12_t_ids;
	std::vector<std::array<int, 2>>& new_edges = buffers.new_edges;

	while( !es_queue.empty() )
	{
		std::array<int, 2> v_ids = es_queue.top().v_ids;
		es_queue.pop();

		if( tet_vertices[ v_ids[ 0 ] ].isFreezed() && tet_vertices[ v_ids[ 1 ] ].isFreezed() )
			continue;

		n12_t_ids.clear();
		setIntersection( tet_vertices[ v_ids[ 0 ] ].connTets, tet_vertices[ v_ids[ 1 ] ].connTets, n12_t_ids );
		if( !is_swappable( v_ids[ 0 ], v_ids[ 1 ], n12_t_ids ) )
			continue;

		while( !es_queue.empty() )
		{
			if( es_queue.top().v_ids == v_ids )
				es_queue.pop();
			else
				break;
		}

		bool is_success = false;
		new_edges.clear();
		if( n12_t_ids.size() == 3 && remove_an_edge_32( mesh, v_ids[ 0 ], v_ids[ 1 ], n12_t_ids, new_edges ) )
		{
			suc_counter3++;
			is_success = true;
			//            output_info(mesh);
		}
		if( n12_t_ids.size() == 4 && remove_an_edge_44( mesh, v_ids[ 0 ], v_ids[ 1 ], n12_t_ids, new_edges ) )
		{
			suc_counter4++;
			is_success = true;
			//            output_info(mesh);
		}
		if( n12_t_ids.size() == 5 && remove_an_edge_56( mesh, v_ids[ 0 ], v_ids[ 1 ], n12_t_ids, new_edges ) )
		{
			suc_counter5++;
			is_success = true;
			//            output_info(mesh);
		}

		for( auto& e : new_edges )
		{
			es_queue.push( ElementInQueue( e, get_edge_length_2( mesh, e[ 0 ], e[ 1 ] ) ) );
		}

		counter++;
	}

	mesh.logger().logDebug( "success3 = %i, success4 = %i, success5 = %i", suc_counter3, suc_counter4, suc_counter5 );
	mesh.logger().logDebug( "success = %i ( %i )", suc_counter3 + suc_counter4 + suc_counter5, counter );
}

namespace
{
	using namespace floatTetWild;

	inline bool isInverted( const std::vector<MeshVertex>& vertices, const Vector4i& tetra )
	{
		return is_inverted( vertices[ tetra[ 0 ] ], vertices[ tetra[ 1 ] ], vertices[ tetra[ 2 ] ], vertices[ tetra[ 3 ] ] );
	}

	inline double getQuality( const std::vector<MeshVertex>& vertices, const Vector4i& tetra )
	{
		return get_quality( vertices[ tetra[ 0 ] ], vertices[ tetra[ 1 ] ], vertices[ tetra[ 2 ] ], vertices[ tetra[ 3 ] ] );
	}
}  // namespace

// Requires 3 elements in old_t_ids, on success produces up to 12 new edges
bool floatTetWild::remove_an_edge_32( Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges )
{
	if( old_t_ids.size() != 3 )
		return false;

	auto& tet_vertices = mesh.tet_vertices;
	auto& tets = mesh.tets;

	////construct
	std::array<int, 2> v_ids;
	std::array<MeshTet, 2> new_tets;
	std::array<int, 2> t_ids;
	int cnt = 0;
	for( int i = 0; i < 4; i++ )
	{
		if( tets[ old_t_ids[ 0 ] ][ i ] != v1_id && tets[ old_t_ids[ 0 ] ][ i ] != v2_id )
			v_ids[ cnt++ ] = tets[ old_t_ids[ 0 ] ][ i ];
	}
	int i = tets[ old_t_ids[ 1 ] ].find( v_ids[ 0 ] );
	if( i >= 0 )
	{
		new_tets[ 0 ] = tets[ old_t_ids[ 1 ] ];
		new_tets[ 1 ] = tets[ old_t_ids[ 2 ] ];
		t_ids = { { old_t_ids[ 1 ], old_t_ids[ 2 ] } };
	}
	else
	{
		new_tets[ 0 ] = tets[ old_t_ids[ 2 ] ];
		new_tets[ 1 ] = tets[ old_t_ids[ 1 ] ];
		t_ids = { { old_t_ids[ 2 ], old_t_ids[ 1 ] } };
	}
	i = new_tets[ 0 ].find( v1_id );
	new_tets[ 0 ][ i ] = v_ids[ 1 ];
	i = new_tets[ 1 ].find( v2_id );
	new_tets[ 1 ][ i ] = v_ids[ 0 ];

	////check
	for( auto& t : new_tets )
		if( isInverted( tet_vertices, t.indices ) )
			return false;

	std::array<double, 2> new_qs;
	double old_max_quality = 0;
	for( int t_id : old_t_ids )
	{
		if( tets[ t_id ].quality > old_max_quality )
			old_max_quality = tets[ t_id ].quality;
	}
	for( size_t i = 0; i < 2; i++ )
	{
		double q = getQuality( tet_vertices, new_tets[ i ].indices );
		if( q >= old_max_quality )	// or use > ???
			return false;
		new_qs[ i ] = q;
	}

	// ==== real update ====
	constexpr size_t maxLength = 6;
	size_t tempItems = 0;
	std::array<std::array<int, 3>, maxLength> fs;
	std::array<int8_t, maxLength> is_sf_fs;
	std::array<int8_t, maxLength> sf_tags;
	std::array<int8_t, maxLength> is_bx_fs;

	for( int i = 0; i < old_t_ids.size(); i++ )
	{
		for( int j = 0; j < 4; j++ )
		{
			if( tets[ old_t_ids[ i ] ][ j ] == v1_id || tets[ old_t_ids[ i ] ][ j ] == v2_id )
			{
				assert( tempItems < maxLength );

				std::array<int, 3>& tmp = fs[ tempItems ];
				const auto& oldTet = tets[ old_t_ids[ i ] ];
				tmp[ 0 ] = oldTet.indices[ mod4( j + 1 ) ];
				tmp[ 1 ] = oldTet.indices[ mod4( j + 2 ) ];
				tmp[ 2 ] = oldTet.indices[ mod4( j + 3 ) ];
				std::sort( tmp.begin(), tmp.end() );

				is_sf_fs[ tempItems ] = oldTet.is_surface_fs[ j ];
				sf_tags[ tempItems ] = oldTet.surface_tags[ j ];
				is_bx_fs[ tempItems ] = oldTet.is_bbox_fs[ j ];

				tempItems++;
			}
		}
	}

	tets[ old_t_ids[ 0 ] ].is_removed = true;
	tets[ t_ids[ 0 ] ] = new_tets[ 0 ];	 // v2
	tets[ t_ids[ 1 ] ] = new_tets[ 1 ];	 // v1

	for( int i = 0; i < 2; i++ )
		tets[ t_ids[ i ] ].quality = new_qs[ i ];

	for( int i = 0; i < 4; i++ )
	{
		if( tets[ t_ids[ 0 ] ][ i ] != v2_id )
		{
			std::array<int, 3> tmp = { { tets[ t_ids[ 0 ] ][ mod4( i + 1 ) ], tets[ t_ids[ 0 ] ][ mod4( i + 2 ) ], tets[ t_ids[ 0 ] ][ mod4( i + 3 ) ] } };
			std::sort( tmp.begin(), tmp.end() );
			auto it = std::find( fs.begin(), fs.begin() + tempItems, tmp );
			tets[ t_ids[ 0 ] ].is_surface_fs[ i ] = is_sf_fs[ it - fs.begin() ];
			tets[ t_ids[ 0 ] ].surface_tags[ i ] = sf_tags[ it - fs.begin() ];
			tets[ t_ids[ 0 ] ].is_bbox_fs[ i ] = is_bx_fs[ it - fs.begin() ];
		}
		else
		{
			tets[ t_ids[ 0 ] ].is_surface_fs[ i ] = NOT_SURFACE;
			tets[ t_ids[ 0 ] ].surface_tags[ i ] = NO_SURFACE_TAG;
			tets[ t_ids[ 0 ] ].is_bbox_fs[ i ] = NOT_BBOX;
		}

		if( tets[ t_ids[ 1 ] ][ i ] != v1_id )
		{
			std::array<int, 3> tmp = { { tets[ t_ids[ 1 ] ][ mod4( i + 1 ) ], tets[ t_ids[ 1 ] ][ mod4( i + 2 ) ], tets[ t_ids[ 1 ] ][ mod4( i + 3 ) ] } };
			std::sort( tmp.begin(), tmp.end() );
			auto it = std::find( fs.begin(), fs.begin() + tempItems, tmp );
			tets[ t_ids[ 1 ] ].is_surface_fs[ i ] = is_sf_fs[ it - fs.begin() ];
			tets[ t_ids[ 1 ] ].surface_tags[ i ] = sf_tags[ it - fs.begin() ];
			tets[ t_ids[ 1 ] ].is_bbox_fs[ i ] = is_bx_fs[ it - fs.begin() ];
		}
		else
		{
			tets[ t_ids[ 1 ] ].is_surface_fs[ i ] = NOT_SURFACE;
			tets[ t_ids[ 1 ] ].surface_tags[ i ] = NO_SURFACE_TAG;
			tets[ t_ids[ 1 ] ].is_bbox_fs[ i ] = NOT_BBOX;
		}
	}

	tet_vertices[ v_ids[ 0 ] ].connTets.remove( old_t_ids[ 0 ] );
	tet_vertices[ v_ids[ 1 ] ].connTets.remove( old_t_ids[ 0 ] );

	tet_vertices[ v_ids[ 0 ] ].connTets.add( t_ids[ 1 ] );
	tet_vertices[ v_ids[ 1 ] ].connTets.add( t_ids[ 0 ] );

	tet_vertices[ v1_id ].connTets.remove( old_t_ids[ 0 ] );
	tet_vertices[ v2_id ].connTets.remove( old_t_ids[ 0 ] );

	tet_vertices[ v1_id ].connTets.remove( t_ids[ 0 ] );
	tet_vertices[ v2_id ].connTets.remove( t_ids[ 1 ] );

	new_edges.reserve( new_tets.size() * 6 );
	for( int i = 0; i < new_tets.size(); i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			std::array<int, 2> e = { { new_tets[ i ][ 0 ], new_tets[ i ][ j + 1 ] } };
			if( e[ 0 ] > e[ 1 ] )
				std::swap( e[ 0 ], e[ 1 ] );
			new_edges.push_back( e );
			e = { { new_tets[ i ][ j + 1 ], new_tets[ i ][ mod3( j + 1 ) + 1 ] } };
			if( e[ 0 ] > e[ 1 ] )
				std::swap( e[ 0 ], e[ 1 ] );
			new_edges.push_back( e );
		}
	}
	vector_unique( new_edges );

	return true;
}

// Requires 4 elements in old_t_ids, on success produces up to 24 new edges
bool floatTetWild::remove_an_edge_44( Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges )
{
	constexpr int N = 4;
	if( old_t_ids.size() != N )
		return false;

	auto& tet_vertices = mesh.tet_vertices;
	auto& tets = mesh.tets;

	////construct
	std::array<std::array<int, 3>, N> n12_es;
	for( int i = 0; i < N; i++ )
	{
		std::array<int, 3>& e = n12_es[ i ];

		int cnt = 0;
		for( int j = 0; j < 4; j++ )
			if( tets[ old_t_ids[ i ] ][ j ] != v1_id && tets[ old_t_ids[ i ] ][ j ] != v2_id )
				e[ cnt++ ] = tets[ old_t_ids[ i ] ][ j ];

		assert( cnt == 2 );
		e[ 2 ] = old_t_ids[ i ];
	}

	SmallBuffer<int, 5> n12_v_ids;
	SmallBuffer<int, 5> n12_t_ids;
	n12_v_ids.push_back( n12_es[ 0 ][ 0 ] );
	n12_v_ids.push_back( n12_es[ 0 ][ 1 ] );
	n12_t_ids.push_back( n12_es[ 0 ][ 2 ] );
	Bitset32<N> is_visited( 1 );
	for( int i = 0; i < N - 2; i++ )	// 2 iterations
	{
		for( int j = 0; j < N; j++ )	// 4 iterations
		{
			if( !is_visited[ j ] )
			{
				if( n12_es[ j ][ 0 ] == n12_v_ids.back() )
				{
					is_visited.set( j );
					n12_v_ids.push_back( n12_es[ j ][ 1 ] );
				}
				else if( n12_es[ j ][ 1 ] == n12_v_ids.back() )
				{  // else if!!!!!!!!!!
					is_visited.set( j );
					n12_v_ids.push_back( n12_es[ j ][ 0 ] );
				}
				if( is_visited[ j ] )
				{
					n12_t_ids.push_back( n12_es[ j ][ 2 ] );
					break;
				}
			}
		}
	}
	n12_t_ids.push_back( n12_es[ is_visited.firstFalseIndex() ][ 2 ] );

	////check
	bool is_valid = false;
	std::array<Vector4i, 4> new_tets;
	std::array<int, 4> tags;
	std::array<int, 2> v_ids;
	std::array<double, 4> new_qs;
	double old_max_quality = 0;
	double new_max_quality = 0;
	for( int t_id : old_t_ids )
		if( tets[ t_id ].quality > old_max_quality )
			old_max_quality = tets[ t_id ].quality;

	for( int i = 0; i < 2; i++ )
	{
		std::array<Vector4i, 4> tmp_new_tets;
		std::array<int, 4> tmp_tags;
		std::array<int, 2> tmp_v_ids;
		tmp_v_ids = { { n12_v_ids[ 0 + i ], n12_v_ids[ 2 + i ] } };
		bool is_break = false;
		for( int j = 0; j < N; j++ )
		{
			auto t = tets[ old_t_ids[ j ] ];
			int ii = t.find( tmp_v_ids[ 0 ] );
			if( ii >= 0 )
			{
				int jt = t.find( v2_id );
				t[ jt ] = tmp_v_ids[ 1 ];
				tmp_tags[ j ] = 1;
			}
			else
			{
				int jt = t.find( v1_id );
				t[ jt ] = tmp_v_ids[ 0 ];
				tmp_tags[ j ] = 0;
			}
			if( isInverted( tet_vertices, t.indices ) )
			{
				is_break = true;
				break;
			}
			tmp_new_tets[ j ] = t.indices;
		}
		if( is_break )
			continue;

		std::array<double, 4> tmp_new_qs;
		for( int j = 0; j < N; j++ )
		{
			double q = getQuality( tet_vertices, tmp_new_tets[ j ] );
			if( q >= old_max_quality )
			{
				is_break = true;
				break;
			}
			if( q > new_max_quality )
				new_max_quality = q;
			tmp_new_qs[ j ] = q;
		}
		if( is_break )
			continue;

		is_valid = true;
		old_max_quality = new_max_quality;
		new_tets = tmp_new_tets;
		tags = tmp_tags;
		new_qs = tmp_new_qs;
		v_ids = tmp_v_ids;
	}
	if( !is_valid )
		return false;

	// ==== real update ====
	constexpr size_t maxLength = 8;
	size_t tempItems = 0;
	std::array<std::array<int, 3>, maxLength> fs;
	std::array<int8_t, maxLength> is_sf_fs;
	std::array<int8_t, maxLength> sf_tags;
	std::array<int8_t, maxLength> is_bx_fs;

	for( int i = 0; i < old_t_ids.size(); i++ )
	{
		const auto& oldTet = tets[ old_t_ids[ i ] ];
		for( int j = 0; j < 4; j++ )
		{
			if( oldTet.indices[ j ] == v1_id || oldTet.indices[ j ] == v2_id )
			{
				std::array<int, 3>& tmp = fs[ tempItems ];
				tmp[ 0 ] = oldTet.indices[ mod4( j + 1 ) ];
				tmp[ 1 ] = oldTet.indices[ mod4( j + 2 ) ];
				tmp[ 2 ] = oldTet.indices[ mod4( j + 3 ) ];
				std::sort( tmp.begin(), tmp.end() );

				is_sf_fs[ tempItems ] = oldTet.is_surface_fs[ j ];
				sf_tags[ tempItems ] = oldTet.surface_tags[ j ];
				is_bx_fs[ tempItems ] = oldTet.is_bbox_fs[ j ];

				tempItems++;
				assert( tempItems <= maxLength );
			}
		}
	}

	for( int j = 0; j < new_tets.size(); j++ )
	{
		if( tags[ j ] == 0 )
		{
			tet_vertices[ v1_id ].connTets.remove( old_t_ids[ j ] );
			tet_vertices[ v_ids[ 0 ] ].connTets.add( old_t_ids[ j ] );
		}
		else
		{
			tet_vertices[ v2_id ].connTets.remove( old_t_ids[ j ] );
			tet_vertices[ v_ids[ 1 ] ].connTets.add( old_t_ids[ j ] );
		}
		tets[ old_t_ids[ j ] ].indices = new_tets[ j ];
		tets[ old_t_ids[ j ] ].quality = new_qs[ j ];
	}

	const auto fsEnd = fs.begin() + tempItems;
	for( int i = 0; i < old_t_ids.size(); i++ )
	{  // old_t_ids contains new tets
		for( int j = 0; j < 4; j++ )
		{
			tets[ old_t_ids[ i ] ].is_surface_fs[ j ] = NOT_SURFACE;
			tets[ old_t_ids[ i ] ].surface_tags[ j ] = NO_SURFACE_TAG;
			tets[ old_t_ids[ i ] ].is_bbox_fs[ j ] = NOT_BBOX;
			if( tets[ old_t_ids[ i ] ][ j ] == v_ids[ 0 ] || tets[ old_t_ids[ i ] ][ j ] == v_ids[ 1 ] )
			{
				std::array<int, 3> tmp = {
				  { tets[ old_t_ids[ i ] ][ mod4( j + 1 ) ], tets[ old_t_ids[ i ] ][ mod4( j + 2 ) ], tets[ old_t_ids[ i ] ][ mod4( j + 3 ) ] } };
				std::sort( tmp.begin(), tmp.end() );
				auto it = std::find( fs.begin(), fsEnd, tmp );
				tets[ old_t_ids[ i ] ].is_surface_fs[ j ] = is_sf_fs[ it - fs.begin() ];
				tets[ old_t_ids[ i ] ].surface_tags[ j ] = sf_tags[ it - fs.begin() ];
				tets[ old_t_ids[ i ] ].is_bbox_fs[ j ] = is_bx_fs[ it - fs.begin() ];
			}
		}
	}

	////re-push
	new_edges.reserve( new_tets.size() * 6 );
	for( int i = 0; i < new_tets.size(); i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			std::array<int, 2> e = { { new_tets[ i ][ 0 ], new_tets[ i ][ j + 1 ] } };
			if( e[ 0 ] > e[ 1 ] )
				std::swap( e[ 0 ], e[ 1 ] );
			new_edges.push_back( e );
			e = { { new_tets[ i ][ j + 1 ], new_tets[ i ][ mod3( j + 1 ) + 1 ] } };
			if( e[ 0 ] > e[ 1 ] )
				std::swap( e[ 0 ], e[ 1 ] );
			new_edges.push_back( e );
		}
	}
	vector_unique( new_edges );

	return true;
}

namespace
{
	// Compute maximum of 6 FP64 values in memory
	inline double computeMaximum( const std::array<double, 6>& arr )
	{
		__m128d v = _mm_loadu_pd( &arr[ 0 ] );
		v = _mm_max_pd( v, _mm_loadu_pd( &arr[ 2 ] ) );
		v = _mm_max_pd( v, _mm_loadu_pd( &arr[ 4 ] ) );
		v = _mm_max_sd( v, _mm_unpackhi_pd( v, v ) );
		return _mm_cvtsd_f64( v );
	}
}

// Requires 5 elements in old_t_ids
bool floatTetWild::remove_an_edge_56( Mesh& mesh, int v1_id, int v2_id, const std::vector<int>& old_t_ids, std::vector<std::array<int, 2>>& new_edges )
{
	if( old_t_ids.size() != 5 )
		return false;

	auto& tet_vertices = mesh.tet_vertices;
	auto& tets = mesh.tets;

	////construct
	std::array<std::array<int, 3>, 5> n12_es;
	for( int i = 0; i < 5; i++ )
	{
		std::array<int, 3>& e = n12_es[ i ];
		int cnt = 0;
		for( int j = 0; j < 4; j++ )
			if( tets[ old_t_ids[ i ] ][ j ] != v1_id && tets[ old_t_ids[ i ] ][ j ] != v2_id )
				e[ cnt++ ] = tets[ old_t_ids[ i ] ][ j ];
		assert( cnt == 2 );
		e[ 2 ] = old_t_ids[ i ];
	}

	SmallBuffer<int, 6> n12_v_ids;
	SmallBuffer<int, 6> n12_t_ids;
	n12_v_ids.push_back( n12_es[ 0 ][ 0 ] );
	n12_v_ids.push_back( n12_es[ 0 ][ 1 ] );
	n12_t_ids.push_back( n12_es[ 0 ][ 2 ] );

	Bitset32<5> is_visited( 1 );
	for( int i = 0; i < 3; i++ )
	{
		for( int j = 0; j < 5; j++ )
		{
			if( !is_visited[ j ] )
			{
				if( n12_es[ j ][ 0 ] == n12_v_ids.back() )
				{
					is_visited.set( j );
					n12_v_ids.push_back( n12_es[ j ][ 1 ] );
				}
				else if( n12_es[ j ][ 1 ] == n12_v_ids.back() )
				{  // else if!!!!!!!!!!
					is_visited.set( j );
					n12_v_ids.push_back( n12_es[ j ][ 0 ] );
				}
				if( is_visited[ j ] )
				{
					n12_t_ids.push_back( n12_es[ j ][ 2 ] );
					break;
				}
			}
		}
	}
	n12_t_ids.push_back( n12_es[ is_visited.firstFalseIndex() ][ 2 ] );

	////check
	double old_max_quality = 0;
	double new_max_quality = 0;
	for( int t_id : old_t_ids )
		old_max_quality = std::max( old_max_quality, tets[ t_id ].quality );

	std::array<std::array<Scalar, 2>, 10> tet_qs;
	std::array<std::array<Vector4i, 2>, 10> new_tets;

	Bitset32<5> is_v_valid( 0b11111 );
	for( int i = 0; i < n12_v_ids.size(); i++ )
	{
		const int iNext = ( i + 1 ) % 5;
		const int iPrev = ( i - 1 + 5 ) % 5;
		const uint32_t nextPrevMask = ( 1u << iNext ) | ( 1u << iPrev );

		// Old version: if( !is_v_valid[ ( i + 1 ) % 5 ] && !is_v_valid[ ( i - 1 + 5 ) % 5 ] )
		if( !is_v_valid.hasAnyBit( nextPrevMask ) )
			continue;

		std::array<Vector4i, 2> new_ts;
		auto t = tets[ n12_t_ids[ i ] ];
		int it = t.find( v1_id );
		t[ it ] = n12_v_ids[ iPrev ];
		if( isInverted( tet_vertices, t.indices ) )
		{
			is_v_valid.clearBits( nextPrevMask );
			continue;
		}
		new_ts[ 0 ] = t.indices;

		t = tets[ n12_t_ids[ i ] ];
		it = t.find( v2_id );
		t[ it ] = n12_v_ids[ iPrev ];
		if( isInverted( tet_vertices, t.indices ) )
		{
			is_v_valid.clearBits( nextPrevMask );
			continue;
		}
		new_ts[ 1 ] = t.indices;

		new_tets[ i ] = new_ts;
		std::array<double, 2>& qs = tet_qs[ i ];

		for( size_t j = 0; j < 2; j++ )
			qs[ j ] = getQuality( tet_vertices, new_ts[ j ] );
	}
	if( is_v_valid.empty() )
		return false;

	int selected_id = -1;
	for( int i = 0; i < 5; i++ )
	{
		if( !is_v_valid[ i ] )
			continue;

		std::array<Vector4i, 2> new_ts;
		auto t = tets[ n12_t_ids[ ( i + 2 ) % 5 ] ];
		int it = t.find( v1_id );
		t[ it ] = n12_v_ids[ i ];
		if( isInverted( tet_vertices, t.indices ) )
			continue;
		new_ts[ 0 ] = t.indices;
		t = tets[ n12_t_ids[ ( i + 2 ) % 5 ] ];
		it = t.find( v2_id );
		t[ it ] = n12_v_ids[ i ];
		if( isInverted( tet_vertices, t.indices ) )
			continue;
		new_ts[ 1 ] = t.indices;

		std::array<double, 6> qs;
		qs[ 0 ] = getQuality( tet_vertices, new_ts[ 0 ] );
		qs[ 1 ] = getQuality( tet_vertices, new_ts[ 1 ] );

		{
			const std::array<double, 2>& a0 = tet_qs[ ( i + 1 ) % 5 ];
			const std::array<double, 2>& a1 = tet_qs[ ( i - 1 + 5 ) % 5 ];
			qs[ 2 ] = a0[ 0 ];
			qs[ 3 ] = a1[ 0 ];
			qs[ 4 ] = a0[ 1 ];
			qs[ 5 ] = a1[ 1 ];
		}
		new_max_quality = std::max( new_max_quality, computeMaximum( qs ) );
		if( new_max_quality >= old_max_quality )
			continue;

		old_max_quality = new_max_quality;
		selected_id = i;

		tet_qs[ i + 5 ] = std::array<Scalar, 2>( { { qs[ 0 ], qs[ 1 ] } } );
		new_tets[ i + 5 ] = new_ts;
	}
	if( selected_id < 0 )
		return false;

	// ==== real update ====
	constexpr size_t maxLength = 10;
	size_t tempItems = 0;
	std::array<std::array<int, 3>, maxLength> fs;
	std::array<int8_t, maxLength> is_sf_fs;
	std::array<int8_t, maxLength> sf_tags;
	std::array<int8_t, maxLength> is_bx_fs;

	// ---- update on surface #1 ----
	for( int i = 0; i < old_t_ids.size(); i++ )
	{
		const MeshTet& oldTet = tets[ old_t_ids[ i ] ];
		for( int j = 0; j < 4; j++ )
		{
			if( oldTet.indices[ j ] == v1_id || oldTet.indices[ j ] == v2_id )
			{
				std::array<int, 3>& tmp = fs[ tempItems ];
				tmp[ 0 ] = oldTet.indices[ mod4( j + 1 ) ];
				tmp[ 1 ] = oldTet.indices[ mod4( j + 2 ) ];
				tmp[ 2 ] = oldTet.indices[ mod4( j + 3 ) ];
				std::sort( tmp.begin(), tmp.end() );

				is_sf_fs[ tempItems ] = oldTet.is_surface_fs[ j ];
				sf_tags[ tempItems ] = oldTet.surface_tags[ j ];
				is_bx_fs[ tempItems ] = oldTet.is_bbox_fs[ j ];

				tempItems++;
				assert( tempItems <= maxLength );
			}
		}
	}

	std::vector<int> new_t_ids = old_t_ids;
	get_new_tet_slots( mesh, 1, new_t_ids );
	tets[ new_t_ids.back() ].reset();

	for( int i = 0; i < 2; i++ )
	{
		tets[ new_t_ids[ i ] ] = new_tets[ ( selected_id + 1 ) % 5 ][ i ];
		tets[ new_t_ids[ i + 2 ] ] = new_tets[ ( selected_id - 1 + 5 ) % 5 ][ i ];
		tets[ new_t_ids[ i + 4 ] ] = new_tets[ selected_id + 5 ][ i ];

		tets[ new_t_ids[ i ] ].quality = tet_qs[ ( selected_id + 1 ) % 5 ][ i ];
		tets[ new_t_ids[ i + 2 ] ].quality = tet_qs[ ( selected_id - 1 + 5 ) % 5 ][ i ];
		tets[ new_t_ids[ i + 4 ] ].quality = tet_qs[ selected_id + 5 ][ i ];
	}

	// update on_surface -- 2
	const auto fsEnd = fs.begin() + tempItems;
	for( int i = 0; i < new_t_ids.size(); i++ )
	{
		for( int j = 0; j < 4; j++ )
		{
			tets[ new_t_ids[ i ] ].is_surface_fs[ j ] = NOT_SURFACE;
			tets[ new_t_ids[ i ] ].surface_tags[ j ] = NO_SURFACE_TAG;
			tets[ new_t_ids[ i ] ].is_bbox_fs[ j ] = NOT_BBOX;
			if( tets[ new_t_ids[ i ] ][ j ] != v1_id && tets[ new_t_ids[ i ] ][ j ] != v2_id &&
				tets[ new_t_ids[ i ] ][ j ] != n12_v_ids[ ( selected_id + 1 ) % 5 ] && tets[ new_t_ids[ i ] ][ j ] != n12_v_ids[ ( selected_id - 1 + 5 ) % 5 ] )
			{
				std::array<int, 3> tmp = {
				  { tets[ new_t_ids[ i ] ][ ( j + 1 ) % 4 ], tets[ new_t_ids[ i ] ][ ( j + 2 ) % 4 ], tets[ new_t_ids[ i ] ][ ( j + 3 ) % 4 ] } };
				std::sort( tmp.begin(), tmp.end() );
				auto it = std::find( fs.begin(), fsEnd, tmp );
				if( it != fsEnd )
				{
					tets[ new_t_ids[ i ] ].is_surface_fs[ j ] = is_sf_fs[ it - fs.begin() ];
					tets[ new_t_ids[ i ] ].surface_tags[ j ] = sf_tags[ it - fs.begin() ];
					tets[ new_t_ids[ i ] ].is_bbox_fs[ j ] = is_bx_fs[ it - fs.begin() ];
				}
			}
		}
	}

	// update conn_tets
	for( int i = 0; i < n12_v_ids.size(); i++ )
	{
		tet_vertices[ n12_v_ids[ i ] ].connTets.remove( n12_t_ids[ i ] );
		tet_vertices[ n12_v_ids[ i ] ].connTets.remove( n12_t_ids[ ( i - 1 + 5 ) % 5 ] );
	}
	for( int i = 0; i < n12_t_ids.size(); i++ )
	{
		tet_vertices[ v1_id ].connTets.remove( n12_t_ids[ i ] );
		tet_vertices[ v2_id ].connTets.remove( n12_t_ids[ i ] );
	}

	// add
	for( int i = 0; i < new_t_ids.size(); i++ )
	{
		for( int j = 0; j < 4; j++ )
			tet_vertices[ tets[ new_t_ids[ i ] ][ j ] ].connTets.add( new_t_ids[ i ] );
	}

	////re-push
	new_edges.reserve( new_t_ids.size() * 6 );
	for( int i = 0; i < new_t_ids.size(); i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			std::array<int, 2> e = { { tets[ new_t_ids[ i ] ][ 0 ], tets[ new_t_ids[ i ] ][ j + 1 ] } };
			if( e[ 0 ] > e[ 1 ] )
				std::swap( e[ 0 ], e[ 1 ] );
			new_edges.push_back( e );
			e = { { tets[ new_t_ids[ i ] ][ j + 1 ], tets[ new_t_ids[ i ] ][ mod3( j + 1 ) + 1 ] } };
			if( e[ 0 ] > e[ 1 ] )
				std::swap( e[ 0 ], e[ 1 ] );
			new_edges.push_back( e );
		}
	}
	vector_unique( new_edges );

	return true;
}