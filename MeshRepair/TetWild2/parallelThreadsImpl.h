#pragma once
#include "../MeshRepair/API/eGlobalFlags.h"

namespace MeshRepair
{
	// RAII class to set thread count thread_local variable for the current thread.
	struct SetThreadsCountRaii
	{
		// If the eGlobalFlags::UseOpenMP was set, set that number to omp_get_max_threads()
		SetThreadsCountRaii( eGlobalFlags flags );

		// Set that number to the specified value
		SetThreadsCountRaii( uint32_t threadsCount );

		// Reset that number back to 0
		~SetThreadsCountRaii();

		SetThreadsCountRaii( const SetThreadsCountRaii& ) = delete;
		void operator=( const SetThreadsCountRaii& ) = delete;
	};
}  // namespace MeshRepair