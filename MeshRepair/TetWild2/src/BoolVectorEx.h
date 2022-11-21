#pragma once
#include <stdint.h>
#include <array>
#include <vector>
#include <assert.h>
#include <queue>

namespace floatTetWild
{
	// This specialized collection is functionally equivalent to std::vector<bool>, optimized for very specific use cases.
	// Exposes an API to safely set these bits from multiple threads concurrently, without locks.
	// Exposes an API to move indices of the set bits to std::queue<int>, without iterating over the complete collection. Sorting is automatic.
	// Overall, this thing is only good for long collections with relatively few set bits.
	// For a dense collection of bits, that outer index adds overhead without any profit.
	class BoolVectorEx
	{
		// A block of 512 bits; the size in bytes and alignment are both 64, that's cache line size
		struct alignas( 64 ) InnerBlock
		{
			std::array<uint8_t, 64> arr;
			InnerBlock()
			{
			}
		};
		// This vector contains individual bits
		std::vector<InnerBlock> inner;

		// This vector is much smaller, it contains one bit per InnerBlock, i.e. approximately 512 times smaller.
		std::vector<uint64_t> outer;

		// Size of the container
		size_t length = 0;

		template<bool atomically>
		static inline void setBit( void* pointer, size_t index )
		{
			uint8_t* const rdi = ( (uint8_t*)pointer ) + index / 8;
			const uint8_t bit = (uint8_t)1 << ( index % 8 );
			if constexpr( atomically )
			{
				// If you're porting this to Linux, try __atomic_or_fetch() or inline assembly
				_InterlockedOr8( (char*)rdi, (char)bit );
			}
			else
				*rdi |= bit;
		}

		static inline bool getBit( const void* pointer, size_t index )
		{
			const uint8_t val = ( (const uint8_t*)pointer )[ index / 8 ];
			const uint8_t bit = (uint8_t)1 << ( index % 8 );
			return 0 != ( val & bit );
		}

	  public:
		void initEmpty( size_t len );

		// Set a value at the specified index
		void setAt( size_t i )
		{
			assert( i < length );
			setBit<false>( inner.data(), i );
			setBit<false>( outer.data(), i / 512 );
		}

		// Set a value at the specified index, concurrently from multiple threads at once
		void setAtomic( size_t i )
		{
			assert( i < length );
			setBit<true>( inner.data(), i );
			setBit<true>( outer.data(), i / 512 );
		}

		bool operator[]( size_t i ) const
		{
			assert( i < length );
			return getBit( inner.data(), i );
		}

		// Enumerate indices of the set values, push these integers to the queue.
		// The queue will then contain unique and sorted numbers.
		void enqueueSet( std::queue<int>& queue ) const;
	};
}  // namespace floatTetWild