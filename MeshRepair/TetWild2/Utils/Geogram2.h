#pragma once
#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/predicates.h>
#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/mesh/mesh_partition.h>
#include <geogram/points/nn_search.h>
#include "Mesh/TriMeshWrapper.h"
#include "Mesh/TriangleMesh.h"
#include "Mesh/TetrahedralMesh.h"
#include "Mesh/TetraMeshWrapper.h"

namespace GEO2
{
	using GEO::Box;
	using GEO::coord_index_t;
	using GEO::index_t;

	using GEO::geo_sgn;
	using GEO::geo_sqr;

	using GEO::NO_FACET;
	using GEO::vec3;
	using GEO::vector;
	namespace PCK
	{
		using GEO::PCK::orient_3d;
	}
	using GEO::bbox_diagonal;
	using GEO::Delaunay;
	using GEO::Delaunay_var;
	using GEO::MESH_ORDER_MORTON;

	using GEO::cross;
	using GEO::distance;
	using GEO::distance2;
	using GEO::dot;
	using GEO::length;
	using GEO::normalize;

	using GEO::Attribute;
	using GEO::mesh_partition;
	using GEO::MESH_PARTITION_HILBERT;
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
	using TetraMesh = ::TetraMeshWrapper;
	using MeshCellsAABB = ::MeshCellsAABBWrapper;
	// using TetraMesh = ::TetrahedralMesh;

	namespace Geom
	{
		inline const vec3& mesh_vertex( const GEO2::Mesh& M, index_t v )
		{
			return M.getVertex( v );
		}
	}  // namespace Geom

	inline double bbox_diagonal( const GEO2::Mesh& M )
	{
		return M.boxDiagonal();
	}
}  // namespace GEO2