#pragma once
#include <stdint.h>

namespace MeshRepair
{
	// It's too hard to pass that integer in parameters to every place which might use OpenMP
	// The implementation uses a thread_local variable for that integer.
	uint32_t getThreadsCount();

	inline bool shouldUseOpenMP()
	{
		return getThreadsCount() > 1;
	}
}  // namespace MeshRepair