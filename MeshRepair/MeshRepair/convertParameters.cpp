#include "stdafx.h"
#include "convertParameters.h"
#include <omp.h>

namespace
{
	inline bool flag( uint32_t flags, MeshRepair::eRepairFlags bit )
	{
		return 0 != ( flags & bit );
	}
}  // namespace

HRESULT convertParameters( floatTetWild::Parameters& rdi, const MeshRepair::Parameters& rsi )
{
	rdi.ideal_edge_length = rsi.idealEdgeLengthRel;
	rdi.eps_rel = rsi.epsilonRel;
	rdi.stop_energy = rsi.stopEnergy;

	using namespace MeshRepair;
	const uint32_t f = rsi.flags;
	rdi.smooth_open_boundary = flag( f, eRepairFlags::SmoothOpenBoundary );
	rdi.manifold_surface = flag( f, eRepairFlags::ManifoldSurface );
	rdi.coarsen = flag( f, eRepairFlags::Coarsen );
	rdi.disable_filtering = flag( f, eRepairFlags::DisableFiltering );
	rdi.use_floodfill = flag( f, eRepairFlags::UseFloodFill );
	rdi.use_general_wn = flag( f, eRepairFlags::UseGeneralWN );
	rdi.use_input_for_wn = flag( f, eRepairFlags::UseInputForWN );

	rdi.num_threads = 0;
	if( flag( f, eRepairFlags::UseOpenMP ) ) 
	{
		int threads = omp_get_max_threads();
		if( threads > 1 )
			rdi.num_threads = (uint32_t)threads;
	}

	return S_OK;
}