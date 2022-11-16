#pragma once
#include <vector>

namespace floatTetWild
{
	// Global variables from TriangleInsertion.cpp
	struct TriangleInsertionVars
	{
		std::vector<std::array<int, 3>> covered_tet_fs;

		double time_find_cutting_tets = 0;
		double time_find_cutting_tets1 = 0;
		double time_find_cutting_tets2 = 0;
		double time_find_cutting_tets3 = 0;
		double time_find_cutting_tets4 = 0;
		double time_cut_mesh = 0;
		double time_cut_mesh1 = 0;
		double time_cut_mesh2 = 0;
		double time_get_intersecting_edges_and_points = 0;
		double time_subdivide_tets = 0;
		double time_push_new_tets = 0;
		double time_push_new_tets1 = 0;
		double time_push_new_tets2 = 0;
		double time_push_new_tets3 = 0;
		double time_simplify_subdivision_result = 0;
		int cnt_snapped = 0;

		double old_time_find_cutting_tets = 0;
		double old_time_cut_mesh = 0;
		double old_time_get_intersecting_edges_and_points = 0;
		double old_time_subdivide_tets = 0;
		double old_time_push_new_tets = 0;
		double old_time_simplify_subdivision_result = 0;
	};

	struct GlobalVariables
	{
		TriangleInsertionVars triangleInsertion;
	};
}