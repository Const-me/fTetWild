// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once
// To set the parameters related
#include "Types.hpp"
#include <array>
#include <vector>
#include "../Utils/Logger.h"

namespace floatTetWild
{
	class Parameters
	{
	  public:
		Parameters( const MeshRepair::sLoggerSetup& sls )
			: logger( sls )
		{
		}

		bool correct_surface_orientation = false;
		bool is_quiet = false;
		int log_level = 3;	// 2;

		bool smooth_open_boundary = false;
		bool manifold_surface = false;
		bool disable_filtering = false;
		bool use_floodfill = false;
		bool use_general_wn = false;
		bool use_input_for_wn = false;
		bool coarsen = false;

#ifdef NEW_ENVELOPE
		std::vector<double> input_epsr_tags;  // same length as the list of input faces
#endif

		// it decides the scale of the box, presents the deviation of the box from the model
		//( in % of  max((xmax-xmin), (ymax-ymin), (zmax-zmin)) of the input points)
		Scalar box_scale = 1 / 15.0;

		// epsilon presents the tolerance permitted (in % of the box diagonal)
		Scalar eps_rel = 1e-3;

		// initial target edge length at every vertex(in % of the box diagonal)
		Scalar ideal_edge_length_rel = 1 / 20.0;
		Scalar min_edge_len_rel = -1;

		int max_its = 80;
		Scalar stop_energy = 10;
		int stage = 2;

		unsigned int num_threads = std::numeric_limits<unsigned int>::max();

		// The documentation of the command-line parameter says "for debugging usage only"
		// Making it constexpr allows linker to drop large chunk of code in get_max_p() function
		static constexpr int stop_p = -1;

		Vector3 bbox_min;
		Vector3 bbox_max;
		Scalar bbox_diag_length;
		Scalar ideal_edge_length;
		Scalar ideal_edge_length_2;
		Scalar eps_input;
		// The following two numbers are loaded as __m128d vector, please keep them together in that order
		double eps_2;
		double eps;

		Scalar eps_delta;
		Scalar dd;
		Scalar min_edge_length;

		Scalar split_threshold;
		Scalar collapse_threshold;
		Scalar split_threshold_2;
		Scalar collapse_threshold_2;

		Scalar eps_coplanar;
		Scalar eps_2_coplanar;
		// The following two numbers are loaded as __m128d vector, please keep them together in that order
		double eps_2_simplification;
		double eps_simplification;
		Scalar dd_simplification;

		void init( Scalar bbox_diag_l )
		{
			if( stage > 5 )
				stage = 5;

			bbox_diag_length = bbox_diag_l;

			ideal_edge_length = bbox_diag_length * ideal_edge_length_rel;
			ideal_edge_length_2 = ideal_edge_length * ideal_edge_length;

			eps_input = bbox_diag_length * eps_rel;
			dd = eps_input;	 // / stage;
			dd /= 1.5;

			double eps_usable = eps_input - dd / std::sqrt( 3 );
			eps_delta = eps_usable * 0.1;
			eps = eps_usable - eps_delta * ( stage - 1 );

			eps_2 = eps * eps;

			eps_coplanar = eps * 0.2;  // better to set it as eps-related
			if( eps_coplanar > bbox_diag_length * 1e-6 )
				eps_coplanar = bbox_diag_length * 1e-6;
			eps_2_coplanar = eps_coplanar * eps_coplanar;

			eps_simplification = eps * 0.8;
			eps_2_simplification = eps_simplification * eps_simplification;
			dd_simplification = dd / eps * eps_simplification;

			if( min_edge_len_rel < 0 )
				min_edge_len_rel = eps_rel;
			min_edge_length = bbox_diag_length * min_edge_len_rel;

			split_threshold = ideal_edge_length * ( 4 / 3.0 );
			collapse_threshold = ideal_edge_length * ( 4 / 5.0 );
			split_threshold_2 = split_threshold * split_threshold;
			collapse_threshold_2 = collapse_threshold * collapse_threshold;
		}

		const Logger logger;
	};
}  // namespace floatTetWild
