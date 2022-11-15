#pragma once
#include "ElementInQueue.h"
#include "FacetRecursionStack.h"
#include "EdgesSet.h"
#include <queue>

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
	struct CollapseEdgeBuffers
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
		std::vector<bool> is_visited;
		std::queue<int> queue_t_ids;
		std::vector<int> n_t_ids;
	};
}