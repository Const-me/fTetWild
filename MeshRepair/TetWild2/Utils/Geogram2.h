#pragma once
#include "Mesh/TriangleMesh.h"
#include "GeometricPrimitives.h"
#include "NearestSearch.h"
#include "orient3D.h"

namespace GEO2
{
	using coord_index_t = uint8_t;
	using index_t = uint32_t;
	using Mesh = TriangleMesh;

	constexpr uint32_t NO_FACET = ~( (uint32_t)0 );
}  // namespace GEO2