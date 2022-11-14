#pragma once
#include <stdint.h>
#include <vector>

namespace floatTetWild
{
	struct FacetRecursionFrame
	{
		double d;
		uint32_t n, b, e;
	};

	struct alignas( 64 ) FacetRecursionStack
	{
		std::vector<FacetRecursionFrame> stack;
	};

	using FacetRecursionStacks = std::vector<FacetRecursionStack>;
}  // namespace floatTetWild