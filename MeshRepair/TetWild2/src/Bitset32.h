#pragma once
#include <stdint.h>

namespace floatTetWild
{
	// Very similar to std::bitset, too bad std::bitset doesn't supports bit scan
	template<size_t length>
	class Bitset32
	{
		uint32_t bitmap;

	  public:
		Bitset32()
		{
			static_assert( length <= 32 );
			bitmap = 0;
		}

		Bitset32( uint32_t initial )
		{
			assert( 0 == ( initial >> length ) );
			bitmap = initial;
		}

		bool operator[]( size_t pos ) const
		{
			assert( pos < length );
			return 0 != ( bitmap & ( 1u << (uint32_t)pos ) );
		}

		void set( size_t pos )
		{
			assert( pos < length );
			bitmap |= ( 1u << (uint32_t)pos );
		}

		uint32_t firstFalseIndex() const
		{
#ifdef __AVX__
			return _tzcnt_u32( ~bitmap );
#else
#error Without BMI1, you can use _BitScanForward in VC++, or __builtin_ctz in gcc/clang, or std::countr_one if you already use C++/20
#endif	// __AVX__
		}

		void clearBits( uint32_t mask )
		{
			assert( 0 == ( mask >> length ) );
#ifdef __AVX__
			bitmap = _andn_u32( mask, bitmap );
#else
			bitmap &= ~mask;
#endif	// __AVX__
		}

		// True if the bitmap contains all of the bits in the argument
		bool hasAllBits( uint32_t mask ) const
		{
			assert( 0 == ( mask >> length ) );
			return mask == ( bitmap & mask );
		}
		// True if the bitmap contains at least one of the bits in the argument
		bool hasAnyBit( uint32_t mask ) const
		{
			assert( 0 == ( mask >> length ) );
			return 0 != ( bitmap & mask );
		}
		// True if the bitmap is completely empty, not even a single bit
		bool empty() const
		{
			return 0 == bitmap;
		}
	};
}  // namespace floatTetWild