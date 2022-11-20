#pragma once
#include "ElementInQueue.h"
#include "FacetRecursionStack.h"
#include "EdgesSet.h"
#include <queue>
#include "../Utils/Geogram2.h"
#include <map>
#include "TrackedSurfaceChanges.h"
#include "FlatIntMap.h"

namespace floatTetWild
{
	// Temporary data used in find_new_pos function.
	// By the way, 64 is cache line size in bytes.
	// We want different instances of these structures in different cache lines. Otherwise, the performance gonna be ruined by cache coherency protocol.
	struct alignas( 64 ) FindNewPosBuffers
	{
		std::vector<int> js;
		std::vector<std::array<double, 12>> Ts;
		EdgesSet edgesTemp;
	};

	// Temporary buffers used by edge_collapsing_aux function
	struct EdgeCollapsingAuxBuffers
	{
		std::vector<std::array<int, 2>> inf_es;
		std::vector<int> inf_e_tss;
		std::vector<int> tet_tss;
		EdgesSet new_edges;
		std::vector<std::array<int, 2>> tmp_inf_es;
		std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_s> ec_queue;
	};

	// Temporary buffers used by collapse_an_edge function
	struct alignas( 64 ) CollapseEdgeBuffers
	{
		std::vector<int> n12_t_ids;
		std::vector<int> n1_t_ids;
		std::vector<int> js_n1_t_ids;
		std::vector<double> new_qs;
		std::vector<int> n1_v_ids;
		std::vector<int> pair;
	};

	// Temporary buffers used by edge_collapsing function
	struct EdgeCollapsingBuffers
	{
		EdgesSet edges;
	};

	struct FindCuttingTetsBuffers
	{
		std::vector<int> cut_t_ids;
		std::vector<bool> is_visited;
		std::queue<int> queue_t_ids;
		std::vector<int> n_t_ids;
		// Used by findCuttingTetsOmp, to sort triangles produced by different threads
		std::vector<int> queueSortIDs;
	};

	struct IsBoundaryEdgeBuffers
	{
		std::vector<GEO2::vec3> points;
	};

	struct EdgeSwappingBuffers
	{
		// EdgesSet edges;
		std::priority_queue<ElementInQueue, std::vector<ElementInQueue>, cmp_l> es_queue;
		std::vector<int> n12_t_ids;
		std::vector<std::array<int, 2>> new_edges;
	};

	struct SubdivideTetsBuffers
	{
		std::vector<int> n_ids;
		std::vector<std::pair<int, Vector3>> centroids;
		std::vector<std::pair<int, Vector3>> tmp_centroids;
		FlatIntMap map_lv_to_v_id, map_lv_to_c;
	};

	struct InsertOneTriangleBuffers
	{
		std::vector<std::array<int, 3>> covered_tet_fs;
		int cnt_snapped = 0;

		std::vector<Vector3> points;
		FlatEdgeMap map_edge_to_intersecting_point;
		std::vector<int> subdivide_t_ids;
		std::vector<int> tmp;
		std::vector<bool> is_mark_surface;

		std::vector<MeshTet> new_tets;
		TSChanges trackedSurfaceChanges;
		std::vector<int> modified_t_ids;
	};

	struct CutMeshBuffers
	{
		std::vector<int> v_ids;
		// TODO: replace with another container which retains the memory
		std::map<int, int> map_v_ids;
		std::vector<std::array<int, 4>> tets;

		std::vector<double> to_plane_dists;
		std::vector<bool> is_snapped;
		std::vector<bool> is_projected;
	};

	struct alignas( 64 ) InsertionBuffers
	{
		InsertOneTriangleBuffers insertOneTriangle;
		FindCuttingTetsBuffers findCuttingTets;
		CutMeshBuffers cutMesh;
		SubdivideTetsBuffers subdivideTets;
	};
}  // namespace floatTetWild