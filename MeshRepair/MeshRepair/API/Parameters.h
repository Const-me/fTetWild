#pragma once
#include "eGlobalFlags.h"

namespace MeshRepair
{
	enum eRepairFlags : uint32_t
	{
		// Smooth the open boundary
		SmoothOpenBoundary = 1,
		// Force the output to be manifold
		ManifoldSurface = 2,
		// Coarsen the output as much as possible
		Coarsen = 4,
		// Disable filtering out outside elements
		DisableFiltering = 8,
		// Use flood fill to extract interior volume
		UseFloodFill = 0x10,
		// Use general winding number
		UseGeneralWN = 0x20,
		// Use input surface for winding number
		UseInputForWN = 0x40,
		// Skip simplification of the input surface mesh. Makes the result more accurate, at the cost of performance.
		SkipSimplify = 0x80,
		// The idealEdgeLength and epsilon fields are absolute numbers; by default, they're relative to the bounding box diagonal
		LengthsAreAbsolute = 0x100,
		// After making the output mesh, downcast vertices to FP32, and re-index the mesh to remove newly degenerate triangles
		DowncastMeshFp32 = 0x200,
	};

	struct Parameters
	{
		// Ideal length of edges
		// Without LengthsAreAbsolute flag, the value is a multiplier for the diagonal of the bounding box of the input mesh
		// With that flag, the value is absolute number
		double idealEdgeLength = 0.05;

		// Allowed tolerance between input mesh and output mesh
		// Without LengthsAreAbsolute flag, the value is a multiplier for the diagonal of the bounding box of the input mesh
		// With that flag, the value is absolute number
		double epsilon = 1e-3;

		// Stop optimization when max energy is lower than this
		double stopEnergy = 10;

		// See eRepairFlags for the values
		uint32_t flags = eRepairFlags::DowncastMeshFp32;

		bool hasFlag( eRepairFlags f ) const
		{
			return 0 != ( flags & f );
		}
	};
}  // namespace MeshRepair