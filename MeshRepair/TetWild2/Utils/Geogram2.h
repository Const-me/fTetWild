#pragma once
#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/points/nn_search.h>
#include "Mesh/TriMeshWrapper.h"
#include "Mesh/TriangleMesh.h"
#include "GeometricPrimitives.h"
#include "NearestSearch.h"

namespace GEO2
{
	using coord_index_t = uint8_t;
	using index_t = uint32_t;

	using GEO::Delaunay;
	using GEO::Delaunay_var;
	using GEO::vector;
	namespace PCK
	{
		using GEO::PCK::orient_3d;
	}
	constexpr uint32_t NO_FACET = GEO::NO_FACET;

	// using Mesh = TriMeshWrapper;
	using Mesh = TriangleMesh;
}  // namespace GEO2