#pragma once
#include <stdint.h>
#include <vector>

namespace floatTetWild
{
	struct alignas( 16 ) FacetRecursionFrame32
	{
		uint32_t n, b, e;
		float d;

		// Store the values; dVec is assumed to be broadcasted, the method uses W lane, slightly faster to do that way
		void store( __m128i nbe, __m128 dVec )
		{
			__m128 xyz = _mm_castsi128_ps( nbe );
			__m128 res = _mm_blend_ps( xyz, dVec, 0b1000 );
			_mm_store_ps( (float*)&n, res );
		}
	};

	struct alignas( 64 ) FacetRecursionStack
	{
		std::vector<FacetRecursionFrame32> stack32;
	};

	using FacetRecursionStacks = std::vector<FacetRecursionStack>;
}  // namespace floatTetWild