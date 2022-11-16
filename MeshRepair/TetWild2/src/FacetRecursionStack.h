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

		void storeIndices( __m128i nbe )
		{
			_mm_store_si128( (__m128i*)&n, nbe );
		}
	};

	struct alignas( 16 ) FacetRecursionFrame32
	{
		uint32_t n, b, e;
		float d;
	};

	struct alignas( 64 ) FacetRecursionStack
	{
		std::vector<FacetRecursionFrame> stack;
		std::vector<FacetRecursionFrame32> stack32;
	};

	using FacetRecursionStacks = std::vector<FacetRecursionStack>;
}  // namespace floatTetWild