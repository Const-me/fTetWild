// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#pragma once
#include <Eigen/Dense>

namespace floatTetWild
{
	using Scalar = double;
	constexpr double SCALAR_ZERO = 1e-8;
	constexpr double SCALAR_ZERO_2 = 1e-16;
	constexpr double SCALAR_ZERO_3 = 1e-24;

	typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

	typedef Eigen::Matrix<Scalar, 3, 3> Matrix3;

	typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
	typedef Eigen::Matrix<Scalar, 2, 1> Vector2;

	typedef Eigen::Matrix<int, 4, 1> Vector4i;
	typedef Eigen::Matrix<int, 3, 1> Vector3i;
	typedef Eigen::Matrix<int, 2, 1> Vector2i;
}  // namespace floatTetWild

// Insert triangles in parallel, incomplete
#define PARALLEL_TRIANGLES_INSERTION 0