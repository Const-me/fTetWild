#pragma once
#include <stdint.h>
#include <vector>

namespace floatTetWild
{
	struct alignas( 32 ) FacetRecursionFrame
	{
		uint32_t n, b, e;
		uint32_t zzPadding;
		double d;

		inline void storeIndices( __m128i nbe )
		{
			_mm_store_si128( (__m128i*)&n, nbe );
		}
	};

	struct alignas( 64 ) FacetRecursionStack
	{
		std::vector<FacetRecursionFrame> stack;
	};

	using FacetRecursionStacks = std::vector<FacetRecursionStack>;
}  // namespace floatTetWild