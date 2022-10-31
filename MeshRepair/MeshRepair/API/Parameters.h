#pragma once
#include <stdint.h>

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
	};

	struct Parameters
	{
		// Ideal length of edges, relative to the bounding box diagonal
		double idealEdgeLengthRel = 0.05;

		// epsilon presents the tolerance permitted, relative to the bounding box diagonal
		double epsilonRel = 1e-3;

		// Stop optimization when max energy is lower than this
		double stopEnergy = 10;

		// See eRepairFlags for the values
		uint32_t flags = 0;

		// Maximum number of threads used
		uint32_t maxThreads = 0;
	};
}  // namespace MeshRepair