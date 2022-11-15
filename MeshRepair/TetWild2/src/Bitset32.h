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
	};
}  // namespace floatTetWild