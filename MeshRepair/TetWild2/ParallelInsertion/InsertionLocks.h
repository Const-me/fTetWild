#pragma once
#include <mutex>
#include <shared_mutex>

namespace floatTetWild
{
	struct InsertionLocks
	{
		// Locked while querying or changing the size of these 3 shared vectors
		std::mutex mutex;

		// Locked for reading while doing most of the things, locked for writing while reallocating any of the 3 vectors
		std::shared_mutex shared;
	};
}