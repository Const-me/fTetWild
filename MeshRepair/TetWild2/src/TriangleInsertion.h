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
#pragma once
#include "Mesh.h"
#include "AABBWrapper.h"
#include "CutMesh.h"
#include "../external/Rational.h"
#include "../external/Predicates.h"
#include <map>
#include "BoolVector.h"

namespace floatTetWild
{
	using TrackSF = std::vector<std::array<std::vector<int>, 4>>;

	void match_surface_fs( const Mesh& mesh, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, BoolVector& is_face_inserted,
	  TrackSF& track_surface_fs );

	void insert_triangles( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::vector<int>& input_tags, Mesh& mesh,
	  BoolVector& is_face_inserted, AABBWrapper& tree, bool is_again );

	void insert_triangles_aux( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::vector<int>& input_tags,
	  Mesh& mesh, BoolVector& is_face_inserted, AABBWrapper& tree, bool is_again );

	void sort_input_faces(
	  const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const Mesh& mesh, std::vector<int>& sorted_f_ids );

	void push_new_tets( Mesh& mesh, TrackSF& track_surface_fs, std::vector<Vector3>& points, std::vector<MeshTet>& new_tets,
	  TrackSF& new_track_surface_fs, std::vector<int>& modified_t_ids, bool is_again );

	void pushNewTetsParallel( Mesh& mesh, TrackSF& track_surface_fs, std::vector<Vector3>& points, std::vector<MeshTet>& new_tets,
	  TrackSF& new_track_surface_fs, std::vector<int>& modified_t_ids, bool is_again, size_t prevVerts, size_t prevTets );

	void simplify_subdivision_result( int insert_f_id, int input_v_size, Mesh& mesh, AABBWrapper& tree, TrackSF& track_surface_fs );

	bool insert_one_triangle( int f_id, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
	  const std::vector<int>& input_tags, Mesh& mesh, TrackSF& track_surface_fs, AABBWrapper& tree, bool is_again );

	void find_cutting_tets( int f_id, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::array<Vector3, 3>& vs,
	  const Mesh& mesh, std::vector<int>& result, bool is_again, size_t countTets );

	bool subdivide_tets( int insert_f_id, const Mesh& mesh, CutMesh& cut_mesh, std::vector<Vector3>& points,
	  std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point, const TrackSF& track_surface_fs, std::vector<int>& subdivide_t_ids,
	  std::vector<bool>& is_mark_surface, std::vector<MeshTet>& new_tets, TrackSF& new_track_surface_fs, std::vector<int>& modified_t_ids,
	  size_t countVertices );

	void pair_track_surface_fs( const Mesh& mesh, TrackSF& track_surface_fs );

	/// edge
	void find_boundary_edges( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const BoolVector& is_face_inserted,
	  const BoolVector& old_is_face_inserted, std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos, std::vector<bool>& is_on_cut_edges,
	  std::vector<std::array<int, 2>>& b_edges, const Logger& log );

	bool insert_boundary_edges( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
	  std::vector<std::pair<std::array<int, 2>, std::vector<int>>>& b_edge_infos, std::vector<bool>& is_on_cut_edges, TrackSF& track_surface_fs, Mesh& mesh,
	  AABBWrapper& tree, BoolVector& is_face_inserted, bool is_again, std::vector<std::array<int, 3>>& known_surface_fs,
	  std::vector<std::array<int, 3>>& known_not_surface_fs );

	bool insert_boundary_edges_get_intersecting_edges_and_points( const std::vector<std::vector<std::pair<int, int>>>& covered_fs_infos,
	  const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::array<int, 2>& e, const std::vector<int>& n_f_ids,
	  TrackSF& track_surface_fs, Mesh& mesh, std::vector<Vector3>& points, std::map<std::array<int, 2>, int>& map_edge_to_intersecting_point,
	  std::vector<int>& snapped_v_ids, std::vector<std::array<int, 3>>& cut_fs, bool is_again );

	/// other
	bool is_uninserted_face_covered( int uninserted_f_id, const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
	  const std::vector<int>& cut_t_ids, const Mesh& mesh );

	void mark_surface_fs( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const std::vector<int>& input_tags,
	  TrackSF& track_surface_fs, const BoolVector& is_face_inserted, const std::vector<std::array<int, 3>>& known_surface_fs,
	  const std::vector<std::array<int, 3>>& known_not_surface_fs, std::vector<std::array<int, 2>>& b_edges, Mesh& mesh, AABBWrapper& tree );

	int get_opp_t_id( int t_id, int j, const Mesh& mesh );

	void myassert( bool b, const std::string& s );

	Vector3 getNormal( const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, int idx );

	// fortest
	using Vector3_r = Eigen::Matrix<triwild::Rational, 3, 1>;

	eOrientation orient_rational( const Vector3_r& p1, const Vector3_r& p2, const Vector3_r& p3, const Vector3_r& p );
}  // namespace floatTetWild