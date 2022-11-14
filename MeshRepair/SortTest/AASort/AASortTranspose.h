#pragma once
#include <emmintrin.h>

namespace AASort
{
	// Maximum count of 4x4 blocks in the vector which can be transposed with these lookup tables
	// The function returns 64
	size_t maxInplaceBlocks();

	// Transpose the specified count of 4x4 blocks
	// The function assumes each 4x4 block was already transposed with vector shuffles, and only moves complete vectors
	void transposeBlocksInplace( __m128i* const begin, size_t countBlocks );
}