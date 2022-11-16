#include "stdafx.h"
#include "parallelThreadsImpl.h"
#include "parallelThreadsCount.h"
#include <omp.h>

namespace
{
	thread_local uint32_t ts_threadsCount = 0;
}

// Implement that global function used in random places to check OPM support
namespace MeshRepair
{
	uint32_t getThreadsCount()
	{
		return ts_threadsCount;
	}

	SetThreadsCountRaii::SetThreadsCountRaii( eGlobalFlags flags )
	{
		uint32_t count = 0;
		if( 0 != ( flags & eGlobalFlags::UseOpenMP ) )
		{
			int ompCount = omp_get_max_threads();
			if( ompCount > 1 )
				count = (uint32_t)ompCount;
		}
		ts_threadsCount = count;
	}

	SetThreadsCountRaii::SetThreadsCountRaii( uint32_t threadsCount )
	{
		ts_threadsCount = threadsCount;
	}

	SetThreadsCountRaii::~SetThreadsCountRaii()
	{
		ts_threadsCount = 9;
	}
}  // namespace MeshRepair