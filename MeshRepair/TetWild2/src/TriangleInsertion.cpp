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
#include <atomic>
#include "../ParallelInsertion/parallelInsertion.h"

floatTetWild::Vector3 floatTetWild::getNormal( const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, int idx )
{
	const Vector3i& tri = ib[ idx ];
	const double* p0 = vb[ tri[ 0 ] ].data();
	const double* p1 = vb[ tri[ 1 ] ].data();
	const double* p2 = vb[ tri[ 2 ] ].data();

	using namespace AvxMath;
	const __m256d a = loadDouble3( p0 );
	const __m256d b = loadDouble3( p1 );
	const __m256d c = loadDouble3( p2 );

	const __m256d e1 = _mm256_sub_pd( b, c );
	const __m256d e2 = _mm256_sub_pd( a, c );
	const __m256d cp = vector3Cross( e1, e2 );
	const __m256d norm = vector3Normalize( cp );

	Vector3 res;
	storeDouble3( res.data(), norm );
	return res;
}

void floatTetWild::match_surface_fs( const Mesh& mesh, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  BoolVector& is_face_inserted, TrackSF& track_surface_fs )
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
				is_face_inserted.set( f_id );
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
  Mesh& mesh, BoolVector& is_face_inserted, AABBWrapper& tree, bool is_again )
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

	template<class E, class Lambda>
	inline void removeIf( std::vector<E>& v, Lambda lambda )
	{
		v.erase( std::remove_if( v.begin(), v.end(), lambda ), v.end() );
	}

	template<class E>
	bool haveCapacity( const std::vector<E>& vec, size_t extra )
	{
		return vec.size() + extra <= vec.capacity();
	}

	static inline size_t addOneK( size_t i )
	{
		constexpr size_t oneK = 1024;
		i += oneK * 2;		 // add 2k
		i &= ~( oneK - 1 );	 // round down by 1k
		return i;
	}

	static inline size_t growExp( size_t i )
	{
		if( i < 1024 )
			return i;
		// Multiply by approximately 1.618 https://en.wikipedia.org/wiki/Golden_ratio
		i = ( i * 414 ) / 0x100;
		// Round down to 3 most significant bits, mostly for lulz
#ifdef __AVX__
		const size_t lz = _lzcnt_u64( i );
		const size_t mask = 0b111ull << ( 61 - lz );
		return i & mask;
#else
#error Not implemented
#endif
	}

	static size_t newCapacity( size_t oldSize, size_t newSize, size_t oldCap )
	{
		size_t inc = addOneK( newSize );
		size_t e = growExp( newSize );
		return std::max( inc, e );
	}

	template<class E>
	void ensureCapacity( std::vector<E>& vec, size_t extra )
	{
		const size_t oldSize = vec.size();
		const size_t newSize = oldSize + extra;
		const size_t oldCap = vec.capacity();
		const size_t newCap = newCapacity( oldSize, newSize, oldCap );
		assert( newCap > newSize );
		vec.reserve( newCap );
	}
}  // namespace

void floatTetWild::insert_triangles_aux( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  const std::vector<int>& input_tags, Mesh& mesh, BoolVector& is_face_inserted, AABBWrapper& tree, bool is_again )
{
	auto tm = mesh.times.insertTrianglesAux.measure();

	BoolVector old_is_face_inserted = is_face_inserted;	 /// is_face_inserted has been initialized in main

	mesh.logger().logInfo( "triangle insertion start, #f = %zu, #v = %zu, #t = %zu", input_faces.size(), mesh.tet_vertices.size(), mesh.tets.size() );
	/////
	TrackSF track_surface_fs( mesh.tets.size() );
	if( !is_again )
		match_surface_fs( mesh, input_vertices, input_faces, is_face_inserted, track_surface_fs );

	int cnt_matched = is_face_inserted.countTrue();
	mesh.logger().logInfo( "matched #f = %i, uninserted #f = %i", cnt_matched, (int)is_face_inserted.size() - cnt_matched );

	std::vector<int> sorted_f_ids;
	sort_input_faces( input_vertices, input_faces, mesh, sorted_f_ids );
	removeIf( sorted_f_ids, [ & ]( int i ) { return is_face_inserted[ i ]; } );

	std::atomic_int cnt_fail = 0;
	std::atomic_int cnt_total = 0;
	pfnInsertTri pfn = [ & ]( int i )
	{
		assert( !is_face_inserted[ i ] );
		cnt_total++;
		if( insert_one_triangle( i, input_vertices, input_faces, input_tags, mesh, track_surface_fs, tree, is_again ) )
			is_face_inserted.setAtomic( i );
		else
			cnt_fail++;
	};

#if PARALLEL_TRIANGLES_INSERTION
	if( mesh.params.num_threads > 1 )
	{
		ensureCapacity( mesh.tet_vertices, 0 );
		ensureCapacity( mesh.tets, 0 );

		// A good clearance value is 2.0 * maximum element size + a little extra for safety
		// Parallel insertion relies on the fact different threads are modifying different elements
		__m256d ts = mesh.maxTetraSize;
		constexpr double parallelClearance = 2.0 + 1.0 / 16.0;
		__m256d clearance = _mm256_mul_pd( ts, _mm256_set1_pd( parallelClearance ) );
		parallelInsertion( input_vertices, input_faces, sorted_f_ids, clearance, pfn );
	}
	else
	{
		for( int i : sorted_f_ids )
			pfn( i );
	}
#else
	for( int i : sorted_f_ids )
		pfn( i );
#endif

	mesh.logger().logInfo( "insert_one_triangle * n done, #v = %zu, #t = %zu", mesh.tet_vertices.size(), mesh.tets.size() );
	mesh.logger().logInfo( "uninserted #f = %zu/%i", is_face_inserted.countFalse(), (int)is_face_inserted.size() - cnt_matched );

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
	mesh.logger().logInfo( "uninserted #f = %zu/%i", is_face_inserted.countFalse(), (int)is_face_inserted.size() - cnt_matched );

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

	if( is_face_inserted.countFalse() == 0 )
		mesh.is_input_all_inserted = true;
	mesh.logger().logInfo( "#b_edge1 = %zu, #b_edges2 = %zu", b_edges1.size(), b_edges2.size() );
	pausee();
}

namespace
{
	template<class E>
	inline void appendVector( std::vector<E>& rdi, const std::vector<E>& rsi )
	{
		rdi.insert( rdi.end(), rsi.begin(), rsi.end() );
	}

	template<class E>
	inline void sortVector( std::vector<E>& vec )
	{
		std::sort( vec.begin(), vec.end() );
	}
}  // namespace

bool floatTetWild::insert_one_triangle( int insert_f_id, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  const std::vector<int>& input_tags, Mesh& mesh, TrackSF& track_surface_fs, AABBWrapper& tree, bool is_again )
{
	auto tm = mesh.times.insertOneTriangle.measure();
#if PARALLEL_TRIANGLES_INSERTION
	size_t countTets, countVertices;
	{
		std::lock_guard<std::mutex> lock { mesh.locks.mutex };
		countTets = mesh.tets.size();
		countVertices = mesh.tet_vertices.size();
	}
	std::shared_lock lockShared { mesh.locks.shared };
#else
	const size_t countTets = mesh.tets.size();
	const size_t countVertices = mesh.tet_vertices.size();
#endif

	std::array<Vector3, 3> vs = { { input_vertices[ input_faces[ insert_f_id ][ 0 ] ], input_vertices[ input_faces[ insert_f_id ][ 1 ] ],
	  input_vertices[ input_faces[ insert_f_id ][ 2 ] ] } };
	Vector3 n = ( vs[ 1 ] - vs[ 0 ] ).cross( vs[ 2 ] - vs[ 0 ] );
	n.normalize();
	const int t = get_t( vs[ 0 ], vs[ 1 ], vs[ 2 ] );

	auto& insertionBuffers = mesh.insertionBuffers();
	std::vector<int>& cut_t_ids = insertionBuffers.findCuttingTets.cut_t_ids;
	cut_t_ids.clear();
	find_cutting_tets( insert_f_id, input_vertices, input_faces, vs, mesh, cut_t_ids, is_again, countTets );

	if( cut_t_ids.empty() )
		throw std::logic_error( "cut_t_ids.empty()" );

	InsertOneTriangleBuffers& buffers = insertionBuffers.insertOneTriangle;
	CutMesh cut_mesh( mesh, n, vs, insertionBuffers );
	cut_mesh.construct( cut_t_ids );

	if( cut_mesh.snap_to_plane() )
	{
		buffers.cnt_snapped++;
		cut_mesh.project_to_plane( input_vertices.size() );
		cut_mesh.expand_new( cut_t_ids, countTets );
		cut_mesh.project_to_plane( input_vertices.size() );
	}

	std::vector<Vector3>& points = buffers.points;
	points.clear();

	std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point = buffers.map_edge_to_intersecting_point;
	map_edge_to_intersecting_point.clear();

	std::vector<int>& subdivide_t_ids = buffers.subdivide_t_ids;
	subdivide_t_ids.clear();

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
	std::vector<int>& tmp = buffers.tmp;
	tmp.clear();

	std::set_difference( subdivide_t_ids.begin(), subdivide_t_ids.end(), cut_t_ids.begin(), cut_t_ids.end(), std::back_inserter( tmp ) );

	std::vector<bool>& is_mark_surface = buffers.is_mark_surface;
	is_mark_surface.clear();
	is_mark_surface.resize( cut_t_ids.size(), true );

	appendVector( cut_t_ids, tmp );
	is_mark_surface.resize( is_mark_surface.size() + tmp.size(), false );

	std::vector<MeshTet>& new_tets = buffers.new_tets;
	new_tets.clear();

	auto& tracked = buffers.trackedSurfaceChanges;
	tracked.clear();

	std::vector<int>& modified_t_ids = buffers.modified_t_ids;
	modified_t_ids.clear();

	if( !subdivide_tets(
		  insert_f_id, mesh, cut_mesh, points, map_edge_to_intersecting_point, cut_t_ids, is_mark_surface, new_tets, tracked, modified_t_ids, countVertices ) )
	{
		if( is_again )
		{
			if( is_uninserted_face_covered( insert_f_id, input_vertices, input_faces, cut_t_ids, mesh ) )
				return true;
		}
		mesh.logger().logWarning( "FAIL subdivide_tets" );
		return false;
	}
#if PARALLEL_TRIANGLES_INSERTION
	lockShared.unlock();
	pushNewTetsParallel( mesh, track_surface_fs, points, new_tets, tracked, modified_t_ids, is_again, countVertices, countTets );
#else
	push_new_tets( mesh, track_surface_fs, points, new_tets, tracked, modified_t_ids, is_again );
#endif

#if PARALLEL_TRIANGLES_INSERTION
	lockShared.lock();
#endif
	simplify_subdivision_result( insert_f_id, input_vertices.size(), mesh, tree, track_surface_fs );
	return true;
}

namespace
{
	static void applyTrackedInserts( TrackSF& track_surface_fs, const TSChanges& tracked, size_t countModified )
	{
		const size_t newSf = tracked.size() - countModified;
		const size_t oldSf = track_surface_fs.size();
		track_surface_fs.resize( oldSf + newSf );
		for( size_t i = 0; i < newSf; i++ )
			tracked[ i + countModified ].apply( track_surface_fs[ i + oldSf ], track_surface_fs );
	}

	static void applyTrackedChanges( TrackSF& track_surface_fs, const TSChanges& tracked, const std::vector<int>& modified_t_ids )
	{
		// If this code will show in the profiler, possible to optimize into a series of in-place moves
		TrackSF newTrackSf;
		newTrackSf.resize( modified_t_ids.size() );
		for( size_t i = 0; i < modified_t_ids.size(); i++ )
			tracked[ i ].apply( newTrackSf[ i ], track_surface_fs );

		// Move these vectors to the destination
		for( size_t i = 0; i < modified_t_ids.size(); i++ )
		{
			const int id = modified_t_ids[ i ];
			track_surface_fs[ id ] = std::move( newTrackSf[ i ] );
		}
	}

	static void applyTracked( TrackSF& track_surface_fs, const TSChanges& tracked, const std::vector<int>& modified_t_ids )
	{
		// Produce new elements of the tracked surface first, this way they view the unmodified values
		applyTrackedInserts( track_surface_fs, tracked, modified_t_ids.size() );

		// Apply edits to the tracked surface
		applyTrackedChanges( track_surface_fs, tracked, modified_t_ids );
	}
}  // namespace

void floatTetWild::push_new_tets( Mesh& mesh, TrackSF& track_surface_fs, std::vector<Vector3>& points, std::vector<MeshTet>& new_tets, const TSChanges& tracked,
  std::vector<int>& modified_t_ids, bool is_again )
{
	// Copy positions of the new vertices
	const int old_v_size = mesh.tet_vertices.size();
	mesh.tet_vertices.resize( mesh.tet_vertices.size() + points.size() );
	for( int i = 0; i < points.size(); i++ )
		mesh.tet_vertices[ old_v_size + i ].pos = points[ i ];

	// Apply changes to the tracked surface
	applyTracked( track_surface_fs, tracked, modified_t_ids );

	// Change mesh topology for the modified elements
	for( int i = 0; i < modified_t_ids.size(); i++ )
	{
		const int id = modified_t_ids[ i ];
		auto& tet = mesh.tets[ id ];
		for( int j = 0; j < 4; j++ )
			mesh.tet_vertices[ tet.indices[ j ] ].connTets.remove( id );

		tet = new_tets[ i ];

		for( int j = 0; j < 4; j++ )
			mesh.tet_vertices[ tet.indices[ j ] ].connTets.add( id );
	}

	for( int i = modified_t_ids.size(); i < new_tets.size(); i++ )
		for( int j = 0; j < 4; j++ )
			mesh.tet_vertices[ new_tets[ i ][ j ] ].connTets.add( mesh.tets.size() + i - modified_t_ids.size() );

	mesh.tets.insert( mesh.tets.end(), new_tets.begin() + modified_t_ids.size(), new_tets.end() );

	modified_t_ids.clear();
}

#if PARALLEL_TRIANGLES_INSERTION
namespace
{
	void offsetNewTets( std::vector<MeshTet>& tets, size_t prevVerts, size_t newVerts )
	{
		if( prevVerts == newVerts || tets.empty() )
			return;

		assert( prevVerts < newVerts );
		const __m128i off = _mm_set1_epi32( (int)( newVerts - prevVerts ) );
		const __m128i pv = _mm_set1_epi32( (int)prevVerts );
		for( MeshTet& tet : tets )
		{
			__m128i* const pointer = (__m128i*)tet.indices.data();
			__m128i v = _mm_loadu_si128( pointer );
			// Compare integers for i < prevVerts
			// When true, the element references a pre-existing vertex, we need to keep such values
			// When false, the element references a new vertex, need to apply offset to compensate for multiple threads changing shared state in parallel
			__m128i isOldVertex = _mm_cmplt_epi32( v, pv );

			__m128i withOffset = _mm_add_epi32( v, off );
			v = _mm_blendv_epi8( withOffset, v, isOldVertex );
			_mm_storeu_si128( pointer, v );
		}
	}

	void offsetTetIDs( std::vector<int>& vec, size_t prevTets, size_t newTets )
	{
		if( prevTets == newTets || vec.empty() )
			return;
		assert( prevTets < newTets );

		const int off = (int)( newTets - prevTets );
		const int pv = (int)prevTets;

		for( int& ref : vec )
		{
			const int i = ref;
			const int withOffset = i + off;
			ref = ( i < pv ) ? i : withOffset;
		}
	}
}  // namespace

void floatTetWild::pushNewTetsParallel( Mesh& mesh, TrackSF& track_surface_fs, std::vector<Vector3>& points, std::vector<MeshTet>& new_tets,
  const TSChanges& tracked, std::vector<int>& modified_t_ids, bool is_again, size_t prevVerts, size_t prevTets )
{
	std::lock_guard<std::mutex> lock { mesh.locks.mutex };

	// Resize the buffers
	size_t vertsIndex, tetsIndex;
	{
		// Compute count of new things
		const size_t newVerts = points.size();
		const size_t newTets = new_tets.size() - modified_t_ids.size();

		// Ensure the capacity
		const bool capVerts = haveCapacity( mesh.tet_vertices, newVerts );
		const bool capTets = haveCapacity( mesh.tets, newTets );

		const bool haveAllCapacity = capVerts && capTets;
		if( !haveAllCapacity )
		{
			// This line kills the concurrency, all threads are keeping shared locks on that shared_mutex
			std::unique_lock lockUnique( mesh.locks.shared );
			// Make sure it's unlikely to happen again
			ensureCapacity( mesh.tet_vertices, newVerts );
			ensureCapacity( mesh.tets, newTets );
		}

		vertsIndex = mesh.tet_vertices.size();
		tetsIndex = mesh.tets.size();

		// Resize the vectors
		mesh.tet_vertices.resize( vertsIndex + newVerts );
		mesh.tets.resize( tetsIndex + newTets );
	}

	// We might need to apply some offsets, to account for other threads adding other things to the mesh
	offsetNewTets( new_tets, prevVerts, vertsIndex );
	offsetTetIDs( modified_t_ids, prevTets, tetsIndex );

	applyTracked( track_surface_fs, tracked, modified_t_ids );

	std::shared_lock lockShared { mesh.locks.shared };

	// Copy positions of the new vertices
	for( size_t i = 0; i < points.size(); i++ )
		mesh.tet_vertices[ vertsIndex + i ].pos = points[ i ];

	// Change mesh topology for the modified elements
	for( size_t i = 0; i < modified_t_ids.size(); i++ )
	{
		for( int j = 0; j < 4; j++ )
			mesh.tet_vertices[ mesh.tets[ modified_t_ids[ i ] ][ j ] ].connTets.remove( modified_t_ids[ i ] );

		mesh.tets[ modified_t_ids[ i ] ] = new_tets[ i ];
		for( int j = 0; j < 4; j++ )
			mesh.tet_vertices[ mesh.tets[ modified_t_ids[ i ] ][ j ] ].connTets.add( modified_t_ids[ i ] );
	}

	// Change mesh topology for the new elements
	for( size_t i = modified_t_ids.size(); i < new_tets.size(); i++ )
		for( int j = 0; j < 4; j++ )
			mesh.tet_vertices[ new_tets[ i ][ j ] ].connTets.add( tetsIndex + i - modified_t_ids.size() );

	// Copy the new elements
	std::copy( new_tets.begin() + modified_t_ids.size(), new_tets.end(), mesh.tets.begin() + tetsIndex );
}
#endif

#include "EdgeCollapsing.h"

void floatTetWild::simplify_subdivision_result( int insert_f_id, int input_v_size, Mesh& mesh, AABBWrapper& tree, TrackSF& track_surface_fs )
{
	std::vector<std::array<int, 3>>& covered_tet_fs = mesh.insertionBuffers().insertOneTriangle.covered_tet_fs;
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

	static void findCuttingTetsOmp( const Mesh& mesh, __m256d min_f, __m256d max_f, std::queue<int>& queue_t_ids, std::vector<bool>& is_visited,
	  FindCuttingTetsBuffers& buffers, size_t countTets )
	{
		std::vector<int>& queueSortIDs = buffers.queueSortIDs;
		queueSortIDs.clear();

		const int64_t length = (int64_t)countTets;
#pragma omp parallel for
		for( int64_t t_id = 0; t_id < length; t_id++ )
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
				queueSortIDs.push_back( (int)t_id );
			}
		}

		std::sort( queueSortIDs.begin(), queueSortIDs.end() );
		for( int i : queueSortIDs )
		{
			queue_t_ids.push( i );
			is_visited[ i ] = true;
		}
	}

	static void findCuttingTets( const Mesh& mesh, __m256d min_f, __m256d max_f, std::queue<int>& queue_t_ids, std::vector<bool>& is_visited, size_t countTets )
	{
		for( size_t t_id = 0; t_id < countTets; t_id++ )
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
  const std::array<Vector3, 3>& vs, const Mesh& mesh, std::vector<int>& cut_t_ids, bool is_again, size_t countTets )
{
	auto tm = mesh.times.findCuttingTets.measure();

	FindCuttingTetsBuffers& buffers = mesh.insertionBuffers().findCuttingTets;

	std::vector<bool>& is_visited = buffers.is_visited;
	is_visited.clear();
	is_visited.resize( countTets, false );

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
#if PARALLEL_TRIANGLES_INSERTION
		findCuttingTets( mesh, min_f, max_f, queue_t_ids, is_visited, countTets );
#else
		if( mesh.params.num_threads > 1 )
			findCuttingTetsOmp( mesh, min_f, max_f, queue_t_ids, is_visited, buffers, countTets );
		else
			findCuttingTets( mesh, min_f, max_f, queue_t_ids, is_visited, countTets );
#endif
	}

	while( !queue_t_ids.empty() )
	{
		const int t_id = queue_t_ids.front();
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
}  // namespace

bool floatTetWild::subdivide_tets( int insert_f_id, const Mesh& mesh, CutMesh& cut_mesh, std::vector<Vector3>& points,
  std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point, std::vector<int>& subdivide_t_ids, std::vector<bool>& is_mark_surface,
  std::vector<MeshTet>& new_tets, TSChanges& tracked, std::vector<int>& modified_t_ids, size_t countVertices )
{
	auto tm = mesh.times.subdivideTets.measure();
	auto& insBuffers = mesh.insertionBuffers();
	SubdivideTetsBuffers& buffers = insBuffers.subdivideTets;

	std::vector<std::array<int, 3>>& covered_tet_fs = insBuffers.insertOneTriangle.covered_tet_fs;
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
						tracked.emplace_back().makeCopyWithAppend( t_id, (uint8_t)j, insert_f_id );
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
		const int v_size = countVertices;
		const int vp_size = countVertices + points.size();
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
		std::vector<int>& n_ids = buffers.n_ids;
		n_ids.clear();

		auto get_centroid = [ & ]( const CutTable::Vec4Buffer& config, int lv_id, Vector3& c )
		{
			n_ids.clear();
			for( const auto& tet : config )
			{
				SmallBuffer<int, 4> tmp;
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

				const double volume = Predicates::orient_3d_volume( vs[ 0 ], vs[ 1 ], vs[ 2 ], vs[ 3 ] );
				if( cnt == 0 )
					min_q = volume;
				else if( volume < min_q )
					min_q = volume;
				cnt++;
			}

			return min_q;
		};

		int diag_config_id = 0;
		std::vector<std::pair<int, Vector3>>& centroids = buffers.centroids;
		centroids.clear();
		std::vector<std::pair<int, Vector3>>& tmp_centroids = buffers.tmp_centroids;

		if( !my_diags.empty() )
		{
			const auto& all_diags = CutTable::get_diag_confs( config_id );
			int bestId = -1;
			double bestQuality = 0;
			for( int i = 0; i < all_diags.size(); i++ )
			{
				if( !all_diags[ i ].equal( my_diags.data(), my_diags.data() + my_diags.size() ) )
					continue;

				tmp_centroids.clear();
				const double min_q = check_config( i, tmp_centroids );
				// The original code of this function was building some nested vectors, sorted them by quality, and then only used the last element of the
				// vector Instead of that expensive stuff, we can search for the best quality on the fly. Copying small vectors ain't that expensive.
				if( min_q > bestQuality )
				{
					bestQuality = min_q;
					bestId = i;
					centroids = tmp_centroids;
				}
			}

			if( bestQuality < SCALAR_ZERO_3 )
			{
				// if tet quality is too bad
				return false;
			}

			diag_config_id = bestId;
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
			auto& newTet =
			  new_tets.emplace_back( MeshTet( map_lv_to_v_id[ t[ 0 ] ], map_lv_to_v_id[ t[ 1 ] ], map_lv_to_v_id[ t[ 2 ] ], map_lv_to_v_id[ t[ 3 ] ] ) );
			auto& tsf = tracked.emplace_back();
			tsf.makeEmpty( t_id );
			for( int j = 0; j < 4; j++ )
			{
				if( new_is_surface_fs[ i ][ j ] && is_mark_sf )
				{
					tsf.prepend( (uint8_t)j, insert_f_id );
					auto& cff = covered_tet_fs.emplace_back();
					cff[ 0 ] = newTet.indices[ ( j + 1 ) % 4 ];
					cff[ 1 ] = newTet.indices[ ( j + 3 ) % 4 ];
					cff[ 2 ] = newTet.indices[ ( j + 2 ) % 4 ];
				}

				int old_local_f_id = new_local_f_ids[ i ][ j ];
				if( old_local_f_id < 0 )
					continue;
				tsf.setSourceCopy( j, old_local_f_id );
				newTet.is_bbox_fs[ j ] = mesh.tets[ t_id ].is_bbox_fs[ old_local_f_id ];
				newTet.is_surface_fs[ j ] = mesh.tets[ t_id ].is_surface_fs[ old_local_f_id ];
				newTet.surface_tags[ j ] = mesh.tets[ t_id ].surface_tags[ old_local_f_id ];
			}
		}
		modified_t_ids.push_back( t_id );
	}

	return true;
}

void floatTetWild::pair_track_surface_fs( const Mesh& mesh, TrackSF& track_surface_fs )
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
			const int opp_t_id = get_opp_t_id( t_id, j, mesh );
			if( opp_t_id < 0 )
				continue;
			const auto& tet = mesh.tets[ t_id ].indices;
			const int k = get_local_f_id( opp_t_id, tet[ ( j + 1 ) % 4 ], tet[ ( j + 2 ) % 4 ], tet[ ( j + 3 ) % 4 ], mesh );
			is_visited[ opp_t_id ][ k ] = true;
			//
			sortVector( track_surface_fs[ t_id ][ j ] );
			sortVector( track_surface_fs[ opp_t_id ][ k ] );
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
  const BoolVector& is_face_inserted, const BoolVector& old_is_face_inserted, std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos,
  std::vector<bool>& is_on_cut_edges, std::vector<std::array<int, 2>>& b_edges, const Logger& log )
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
			sortInt2( e );
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

			Vector3 n = getNormal( input_vertices, input_faces, f_id );
			int t = get_t( input_vertices, input_faces, f_id );

			bool is_fine = false;
			for( int k = 0; k < n12_f_ids.size(); k++ )
			{
				if( n12_f_ids[ k ] == f_id )
					continue;
				Vector3 n1 = getNormal( input_vertices, input_faces, n12_f_ids[ k ] );
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
  std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos, std::vector<bool>& is_on_cut_edges, TrackSF& track_surface_fs, Mesh& mesh,
  AABBWrapper& tree, BoolVector& is_face_inserted, bool is_again, std::vector<std::array<int, 3>>& known_surface_fs,
  std::vector<std::array<int, 3>>& known_not_surface_fs )
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
				is_face_inserted.reset( f_id );
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
		CutMesh empty_cut_mesh( mesh, Vector3( 0, 0, 0 ), std::array<Vector3, 3>(), mesh.insertionBuffers() );
		//
		std::vector<MeshTet> new_tets;
		TSChanges tracked;
		std::vector<int> modified_t_ids;
		if( !subdivide_tets( -1, mesh, empty_cut_mesh, points, map_edge_to_intersecting_point, cut_t_ids, is_mark_surface, new_tets, tracked, modified_t_ids,
			  mesh.tet_vertices.size() ) )
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
					is_face_inserted.reset( f_id );
				for( auto& f : cut_fs )
					mark_known_surface_fs( f, KNOWN_NOT_SURFACE );
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
		push_new_tets( mesh, track_surface_fs, points, new_tets, tracked, modified_t_ids, is_again );
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
  TrackSF& track_surface_fs, Mesh& mesh, std::vector<Vector3>& points, std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
  std::vector<int>& snapped_v_ids, std::vector<std::array<int, 3>>& cut_fs, bool is_again )
{
	//    igl::Timer timer;

	auto is_cross = []( eOrientation a, eOrientation b )
	{
		if( ( a == eOrientation::Positive && b == eOrientation::Negative ) || ( a == eOrientation::Negative && b == eOrientation::Positive ) )
			return true;
		return false;
	};

	const int t = get_t( input_vertices, input_faces, n_f_ids.front() );
	std::array<Vector2, 2> evs_2d = { { to_2d( input_vertices[ e[ 0 ] ], t ), to_2d( input_vertices[ e[ 1 ] ], t ) } };
	const Vector3 n = getNormal( input_vertices, input_faces, n_f_ids.front() );
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
  TrackSF& track_surface_fs, const BoolVector& is_face_inserted, const std::vector<std::array<int, 3>>& known_surface_fs,
  const std::vector<std::array<int, 3>>& known_not_surface_fs, std::vector<std::array<int, 2>>& b_edges, Mesh& mesh, AABBWrapper& tree )
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

bool floatTetWild::is_uninserted_face_covered( int uninserted_f_id, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  const std::vector<int>& cut_t_ids, const Mesh& mesh )
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