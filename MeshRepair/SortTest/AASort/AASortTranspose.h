#pragma once
#include <emmintrin.h>

namespace AASort
{
	void transposeResultOuter( __m128i* const begin, size_t countBlocks );
}