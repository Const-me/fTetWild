// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once
// to generate delaunay mesh using the 3d model along with points of the box
#include "Mesh.h"
#include "Types.hpp"
#include "AABBWrapper.h"
#include "BoolVector.h"

namespace floatTetWild
{
	struct FloatTetDelaunay
	{
		// Generate Delaunay meshes with the model vertices and bounding box points.
		static void tetrahedralize( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const AABBWrapper& tree, Mesh& mesh,
		  BoolVector& is_face_inserted );
	};
}  // namespace floatTetWild