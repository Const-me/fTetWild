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
			cycle = &s_permuteInner[ start ];
			cycleEnd = &s_permuteInner[ end ];
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
}  // namespace

namespace AASort
{

	void transposeWithLookupTable( __m128i* const begin, size_t countBlocks )
	{
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

	void transposeWithLookupTable_x2( __m128i* const begin, size_t countBlocks )
	{
		assert( countBlocks <= lookupTableSizeBlocks * lookupTableSizeBlocks );

		// First pass, transpose the 4xN pieces where N is found in the lookup table
		const size_t countOuter = ( countBlocks + lookupTableSizeBlocks - 1 ) / lookupTableSizeBlocks;
		for( size_t i = 0; i < countOuter; i++ )
		{
			const size_t i1 = i * lookupTableSizeBlocks;
			size_t i2 = ( i + 1 ) * lookupTableSizeBlocks;
			i2 = std::min( i2, countBlocks );
			size_t len = i2 - i1;
			transposeWithLookupTable( begin + i1 * 4, len );
		}

		__m128i buffer[ lookupTableSizeBlocks ];


		throw std::exception( "Not implemented" );
	}
}  // namespace AASort

namespace AASort
{
	void transposeResultOuter( __m128i* const begin, size_t countBlocks )
	{
		if( countBlocks <= lookupTableSizeBlocks )
		{
			transposeWithLookupTable( begin, countBlocks );
			return;
		}
		if( countBlocks <= lookupTableSizeBlocks * lookupTableSizeBlocks )
		{
			transposeWithLookupTable_x2( begin, countBlocks );
			return;
		}
		throw std::exception( "Not implemented" );
	}
}  // namespace AASort