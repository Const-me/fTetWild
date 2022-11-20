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
#pragma once
#include "Mesh.h"
#include "AABBWrapper.h"
#include <map>

namespace floatTetWild
{
	class CutMesh
	{
	  public:
		std::vector<int> v_ids;
		std::map<int, int> map_v_ids;
		std::vector<std::array<int, 4>> tets;

		std::vector<double> to_plane_dists;
		std::vector<bool> is_snapped;
		std::vector<bool> is_projected;

		Mesh& mesh;
		const Vector3& p_n;
		const std::array<Vector3, 3>& p_vs;
		const Logger& logger;

	  private:
		CutMeshBuffers& buffers;

	  public:
		CutMesh( Mesh& _mesh, const Vector3& _p_n, const std::array<Vector3, 3>& _p_vs, InsertionBuffers& insb )
			: mesh( _mesh )
			, p_n( _p_n )
			, p_vs( _p_vs )
			, logger( mesh.params.logger )
			, buffers( insb.cutMesh )
		{
			v_ids = std::move( buffers.v_ids );
			v_ids.clear();

			map_v_ids = std::move( buffers.map_v_ids );
			map_v_ids.clear();

			tets = std::move( buffers.tets );
			tets.clear();

			to_plane_dists = std::move( buffers.to_plane_dists );
			to_plane_dists.clear();

			is_snapped = std::move( buffers.is_snapped );
			is_snapped.clear();

			is_projected = std::move( buffers.is_projected );
			is_projected.clear();
		}

		~CutMesh()
		{
			buffers.v_ids = std::move( v_ids );
			buffers.map_v_ids = std::move( map_v_ids );
			buffers.tets = std::move( tets );
			buffers.to_plane_dists = std::move( to_plane_dists );
			buffers.is_snapped = std::move( is_snapped );
			buffers.is_projected = std::move( is_projected );
		}

		void construct( const std::vector<int>& cut_t_ids );

		bool snap_to_plane();

		void expand( std::vector<int>& cut_t_ids );
		void expand_new( std::vector<int>& cut_t_ids, size_t countTets );

		int project_to_plane( int input_vertices_size );

		bool get_intersecting_edges_and_points( std::vector<Vector3>& points, FlatEdgeMap& map_edge_to_intersecting_point, std::vector<int>& subdivide_t_ids );

		void revert_totally_snapped_tets( int a, int b );

		inline bool is_v_on_plane( int lv_id )
		{
			if( is_snapped[ lv_id ] || to_plane_dists[ lv_id ] == 0 )
				return true;
			return false;
		}

		inline Scalar get_to_plane_dist( const Vector3& p )
		{
			return p_n.dot( p - p_vs[ 0 ] );
		}

		bool check();
	};
}  // namespace floatTetWild