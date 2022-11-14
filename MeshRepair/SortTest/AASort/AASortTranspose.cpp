#include "AASortTranspose.h"
#include <stdint.h>
#include <array>
#include <stdexcept>
#include <assert.h>

namespace
{
#include "TransposePermutations.inl"

	constexpr size_t lookupTableSizeBlocks = s_permuteOuter.size() - 1;

	// Iterator to read one cycle of the permutation
	class Cycle
	{
		const uint8_t* cycle;
		const uint8_t* cycleEnd;

	  public:
		__forceinline Cycle( size_t start, size_t end )
		{
			assert( end > start );	// The compressor filters out empty cycles
			cycle = &s_permuteInner[ start ];

			// C++ guarantees the pointers "one past the end" are valid.
			// However, we can't use operator[] of the array to make such a pointer, because it fails with assertion in debug builds
			assert( end <= s_permuteInner.size() );
			cycleEnd = s_permuteInner.data() + end;
		}

		__forceinline bool notFinished() const
		{
			return cycle < cycleEnd;
		}

		__forceinline size_t next()
		{
			assert( notFinished() );
			size_t res = *cycle;
			cycle++;
			return res;
		}
	};

	// Iterator to read the collection of cycles in the permutation
	class CyclesRange
	{
		const uint16_t* cycles;
		const uint16_t* cyclesEnd;

	  public:
		__forceinline CyclesRange( size_t countBlocks )
		{
			assert( countBlocks <= lookupTableSizeBlocks );

			const size_t outerBegin = s_permuteOuter[ countBlocks - 1 ];
			const size_t outerEnd = s_permuteOuter[ countBlocks ];

			// Thus range is sometimes empty, that's why the ">=" in the assert.
			// More specifically, when transposing a 4x4 block, the compressed collection returns an empty range.
			assert( outerEnd >= outerBegin );

			cycles = &s_permuteCycles[ outerBegin ];
			cyclesEnd = &s_permuteCycles[ outerEnd ];
		}

		__forceinline bool notFinished() const
		{
			return cycles < cyclesEnd;
		}

		__forceinline Cycle next()
		{
			assert( notFinished() );

			const size_t cycleStartIndex = *cycles;
			cycles++;
			const size_t cycleEndIndex = *cycles;
			return Cycle { cycleStartIndex, cycleEndIndex };
		}
	};

	void transposeWithLookupTable( __m128i* const begin, size_t countBlocks )
	{
		// The source code looks complicated, but VC++ only made 36 instructions from this function
		// 2 of them are NOP instructions to pad the two loops, which do nothing.

		assert( countBlocks <= lookupTableSizeBlocks );
		__m128i* const sourceEnd = begin + ( countBlocks * 4 );

		// Each iteration of the inner loop does 2 loads: one scalar from the index buffer, and one vector from the array being transposed
		// And only 1 store, a vector into the array being transposed
		// Should be quite efficient because modern CPUs can do twice as many loads per cycle, compared to stores
		CyclesRange cycles { countBlocks };
		while( cycles.notFinished() )
		{
			Cycle cycle = cycles.next();

			// The compressor filters out empty cycles, it guarantees at least 2 elements in the cycle
			__m128i* const ptrFirst = begin + cycle.next();
			assert( ptrFirst < sourceEnd );
			__m128i* ptr = ptrFirst;
			__m128i prevVec = _mm_loadu_si128( ptr );

			while( cycle.notFinished() )
			{
				ptr = begin + cycle.next();
				assert( ptr < sourceEnd );
				const __m128i nextVec = _mm_loadu_si128( ptr );
				_mm_storeu_si128( ptr, prevVec );
				prevVec = nextVec;
			}

			// Complete the cycle
			_mm_storeu_si128( ptrFirst, prevVec );
		}
	}
}  // namespace AASort

namespace AASort
{
	size_t maxInplaceBlocks()
	{
		return lookupTableSizeBlocks;
	}

	void transposeBlocksInplace( __m128i* const begin, size_t countBlocks )
	{
		assert( countBlocks <= lookupTableSizeBlocks );
		transposeWithLookupTable( begin, countBlocks );
	}
}  // namespace AASort