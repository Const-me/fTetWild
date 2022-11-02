#pragma once
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/predicates.h>
#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_partition.h>
#include <geogram/points/nn_search.h>
#include <geogram/mesh/mesh_AABB.h>

namespace GEO2
{
	using GEO::Mesh;
	using GEO::vec3;
	using GEO::Box;
	using GEO::index_t;
	using GEO::vector;
	using GEO::NO_FACET;
	namespace PCK
	{
		using GEO::PCK::orient_3d;
	}
	using GEO::bbox_diagonal;
	using GEO::MESH_ORDER_MORTON;
	using GEO::Delaunay;
	using GEO::Delaunay_var;

	using GEO::normalize;
	using GEO::dot;
	using GEO::cross;
	using GEO::length;
	using GEO::distance;

	using GEO::mesh_partition;
	using GEO::MESH_PARTITION_HILBERT;
	using GEO::Attribute;
	using GEO::NearestNeighborSearch;
	using GEO::NearestNeighborSearch_var;
	using GEO::MeshCellsAABB;

	namespace Geom
	{
		using GEO::Geom::point_triangle_squared_distance;
	}
}