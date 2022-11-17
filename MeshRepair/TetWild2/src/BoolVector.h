#pragma once
#include <stdint.h>
#include <vector>
#include <assert.h>

namespace floatTetWild
{
	// Same as std::vector<bool>, uses popcnt to compute count of set bits
	class BoolVector
	{
		std::vector<uint64_t> vec;
		size_t length = 0;

	  public:
		BoolVector() = default;
		BoolVector( const BoolVector& ) = default;
		BoolVector( BoolVector&& ) = default;

		uint64_t countTrue() const;
		uint64_t countFalse() const
		{
			return length - countTrue();
		}

		size_t size() const
		{
			return length;
		}

		void set( size_t idx )
		{
			assert( idx < length );
			const size_t elt = idx / 64;
			const uint64_t bit = 1ull << ( idx % 64 );
			vec[ elt ] |= bit;
		}

		void setAtomic( size_t idx )
		{
			assert( idx < length );
			const size_t elt = idx / 64;
			const uint64_t bit = 1ull << ( idx % 64 );

			int64_t* const rdi = (int64_t*)&vec[ elt ];
			_InterlockedOr64( rdi, (int64_t)bit );
		}

		void reset( size_t idx )
		{
			assert( idx < length );
			const size_t elt = idx / 64;
			const uint64_t bit = 1ull << ( idx % 64 );
			const uint64_t invBit = ~bit;
			vec[ elt ] &= invBit;
		}

		bool operator[]( size_t idx ) const
		{
			assert( idx < length );
			const size_t elt = idx / 64;
			const uint64_t bit = 1ull << ( idx % 64 );
			return 0 != ( vec[ elt ] & bit );
		}

		void resize( size_t len, bool fillValue );
	};
}  // namespace floatTetWild