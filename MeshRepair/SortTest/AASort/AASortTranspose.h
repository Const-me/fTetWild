#pragma once
#include <emmintrin.h>

namespace AASort
{
	// Maximum count of 4x4 blocks in the vector which can be efficiently transposed in-place
	size_t maxInplaceBlocks();

	void transposeBlocksInplace( __m128i* const begin, size_t countBlocks );
}