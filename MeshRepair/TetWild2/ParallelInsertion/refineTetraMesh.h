#pragma once
#include <GeogramDelaunay.h>
#include "../src/MeshVertex.h"

namespace floatTetWild
{
	// Compute maximum element size of the Delaunay tetrahedralization; when too large, insert extra vertices into the mesh, and run Delaunay once again
	__m256d refineTetraMesh( const Vector3& boxMin, const Vector3& boxMax, iDelaunay& delaunay, std::vector<Vector3>& vertices, __m128i voxels );
}  // namespace floatTetWild