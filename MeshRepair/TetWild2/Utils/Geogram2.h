#pragma once
#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/points/nn_search.h>
#include "Mesh/TriMeshWrapper.h"
#include "Mesh/TriangleMesh.h"

namespace GEO2
{
	using GEO::Box;
	using GEO::coord_index_t;
	using GEO::index_t;

	using GEO::geo_sgn;
	using GEO::geo_sqr;

	using GEO::vec3;
	using GEO::vector;
	using GEO::Delaunay;
	using GEO::Delaunay_var;
	namespace PCK
	{
		using GEO::PCK::orient_3d;
	}
	constexpr uint32_t NO_FACET = GEO::NO_FACET;

	using GEO::cross;
	using GEO::distance;
	using GEO::distance2;
	using GEO::dot;
	using GEO::length;
	using GEO::normalize;
	using GEO::NearestNeighborSearch;
	using GEO::NearestNeighborSearch_var;

	namespace Geom
	{
		using GEO::Geom::distance2;
		using GEO::Geom::point_triangle_squared_distance;
		using GEO::Geom::tetra_signed_volume;
	}  // namespace Geom

	// using Mesh = TriMeshWrapper;
	using Mesh = TriangleMesh;
}  // namespace GEO2