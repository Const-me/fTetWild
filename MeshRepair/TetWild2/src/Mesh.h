// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License v. 2.0. If a copy of the MPL was not distributed with this file, you can obtain
// one at http://mozilla.org/MPL/2.0/
#pragma once
#include "Parameters.h"
#include "Types.hpp"
#include <vector>
#include <array>
#include <unordered_set>
#include <cassert>
#include <queue>
#include "EdgesSet.h"
#include "Random.hpp"
#include "MeshVertex.h"
#include "MeshTet.h"

namespace floatTetWild
{
	constexpr double MAX_ENERGY = 1e50;

	class Mesh
	{
	  public:
		std::vector<MeshVertex> tet_vertices;
		std::vector<MeshTet> tets;
		Parameters params;
		Mesh( const MeshRepair::sLoggerSetup& sls )
			: params( sls )
		{
		}
		const Logger& logger() const
		{
			return params.logger;
		}

		int t_empty_start = 0;
		int v_empty_start = 0;
		bool is_limit_length = true;
		bool is_closed = true;

		bool is_input_all_inserted = false;
		bool is_coarsening = false;

		void one_ring_vertex_coloring( std::vector<Scalar>& colors ) const;

		void one_ring_vertex_sets( const int threshold, std::vector<std::vector<int>>& concurrent_sets, std::vector<int>& serial_set ) const;

		static void one_ring_edge_set( const std::vector<std::array<int, 2>>& edges, const std::vector<bool>& v_is_removed,
		  const std::vector<bool>& f_is_removed, const std::vector<std::unordered_set<int>>& conn_fs, const std::vector<Vector3>& input_vertices,
		  std::vector<int>& safe_set );

		inline int t_empty_size() const
		{
			int cnt = 0;
			for( const auto& t : tets )
			{
				if( t.is_removed )
					cnt++;
			}
			return cnt;
		}

		inline int v_empty_size() const
		{
			int cnt = 0;
			for( const auto& v : tet_vertices )
			{
				if( v.isRemoved() )
					cnt++;
			}
			return cnt;
		}

		inline void reset_t_empty_start()
		{
			t_empty_start = tets.size();
			for( int i = 0; i < tets.size(); i++ )
			{
				if( tets[ i ].is_removed )
				{
					t_empty_start = i;
					break;
				}
			}
		}

		inline void reset_v_empty_start()
		{
			v_empty_start = tet_vertices.size();
			for( int i = 0; i < tet_vertices.size(); i++ )
			{
				if( tet_vertices[ i ].isRemoved() )
				{
					v_empty_start = i;
					break;
				}
			}
		}

		inline int get_v_num() const
		{
			int cnt = 0;
			for( int i = 0; i < tet_vertices.size(); i++ )
			{
				if( !tet_vertices[ i ].isRemoved() )
					cnt++;
			}
			return cnt;
		}

		inline int get_t_num() const
		{
			int cnt = 0;
			for( int i = 0; i < tets.size(); i++ )
			{
				if( !tets[ i ].is_removed )
					cnt++;
			}
			return cnt;
		}

		inline Scalar get_max_energy() const
		{
			Scalar max_energy = 0;
			int cnt = 0;
			for( auto& t : tets )
			{
				if( t.is_removed )
					continue;
				if( t.quality > max_energy )
					max_energy = t.quality;
			}
			return max_energy;
		}

		inline Scalar get_avg_energy() const
		{
			Scalar avg_energy = 0;
			int cnt = 0;
			for( auto& t : tets )
			{
				if( t.is_removed )
					continue;
				avg_energy += t.quality;
				cnt++;
			}
			avg_energy /= cnt;
			return avg_energy;
		}

		// Temporary data used in find_new_pos function.
		// By the way, 64 is cache line size in bytes.
		// We want different instances of these structures in different cache lines. Otherwise, the performance gonna be ruined by cache coherency protocol.
		struct alignas( 64 ) FindNewPosBuffers
		{
			std::vector<int> js;
			std::vector<std::array<double, 12>> Ts;
			EdgesSet edgesTemp;
		};
		std::vector<FindNewPosBuffers> findNewPosBuffers;

		// Temporary buffers used by edge_collapsing_aux function
		struct EdgeCollapsingAuxBuffers
		{
			std::vector<std::array<int, 2>> inf_es;
			std::vector<int> inf_e_tss;
			std::vector<int> tet_tss;
			EdgesSet new_edges;
			std::vector<std::array<int, 2>> tmp_inf_es;
		};
		EdgeCollapsingAuxBuffers edgeCollapsingAuxBuffers;

		// Temporary buffers used by collapse_an_edge function
		struct CollapseEdgeBuffers
		{
			std::vector<int> n12_t_ids;
			std::vector<int> n1_t_ids;
			std::vector<int> js_n1_t_ids;
			std::vector<Scalar> new_qs;
			std::vector<int> n1_v_ids;
			std::vector<int> pair;
		};
		CollapseEdgeBuffers collapseEdgeBuffers;

		// Temporary buffers used by edge_collapsing function
		struct EdgeCollapsingBuffers
		{
			EdgesSet edges;
		};
		EdgeCollapsingBuffers edgeCollapsingBuffers;

		struct FindCuttingTetsBuffers
		{
			std::vector<bool> is_visited;
			std::queue<int> queue_t_ids;
		};
		FindCuttingTetsBuffers findCuttingTetsBuffers;

		void createThreadLocalBuffers();
	};
}  // namespace floatTetWild