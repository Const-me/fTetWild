#pragma once
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <emmintrin.h>
#include <smmintrin.h>
#ifndef _MSC_VER
#define __forceinline inline __attribute__( ( always_inline ) )
#endif

namespace
{
	struct UintTraits
	{
		static __forceinline __m128i min( __m128i a, __m128i b )
		{
			return _mm_min_epu32( a, b );
		}
		static __forceinline __m128i max( __m128i a, __m128i b )
		{
			return _mm_max_epu32( a, b );
		}
		static constexpr uint32_t minValue = 0;
		static constexpr uint32_t maxValue = ~minValue;
	};

	struct IntTraits
	{
		static __forceinline __m128i min( __m128i a, __m128i b )
		{
			return _mm_min_epi32( a, b );
		}
		static __forceinline __m128i max( __m128i a, __m128i b )
		{
			return _mm_max_epi32( a, b );
		}
		static constexpr int minValue = INT_MIN;
		static constexpr int maxValue = INT_MAX;
	};

	// Transpose 4x4 matrix of int32 elements, in 4 vectors
	__forceinline void transpose4x4( __m128i& v0, __m128i& v1, __m128i& v2, __m128i& v3 )
	{
		// https://randombit.net/bitbashing/posts/integer_matrix_transpose_in_sse2.html
		__m128i t0 = _mm_unpacklo_epi32( v0, v1 );
		__m128i t1 = _mm_unpacklo_epi32( v2, v3 );
		__m128i t2 = _mm_unpackhi_epi32( v0, v1 );
		__m128i t3 = _mm_unpackhi_epi32( v2, v3 );

		v0 = _mm_unpacklo_epi64( t0, t1 );
		v1 = _mm_unpackhi_epi64( t0, t1 );
		v2 = _mm_unpacklo_epi64( t2, t3 );
		v3 = _mm_unpackhi_epi64( t2, t3 );
	}

	template<class Traits>
	struct Primitives : public Traits
	{
		// Vertically sort 4 integer vectors
		static __forceinline void sort4( __m128i& a, __m128i& b, __m128i& c, __m128i& d )
		{
			__m128i i01 = Traits::min( a, b );
			__m128i ax01 = Traits::max( a, b );

			__m128i i23 = Traits::min( c, d );
			__m128i ax23 = Traits::max( c, d );

			a = Traits::min( i01, i23 );
			d = Traits::max( ax01, ax23 );

			__m128i m0 = Traits::max( i01, i23 );
			__m128i m1 = Traits::min( ax01, ax23 );
			b = Traits::min( m0, m1 );
			c = Traits::max( m0, m1 );
		}

		// Load 4 vectors, vertically sort them, and store in transposed order
		static __forceinline void copySortedBlock( __m128i* rdi, const void* rsi )
		{
			const __m128i* const src = (const __m128i*)rsi;
			__m128i v0 = _mm_loadu_si128( src );
			__m128i v1 = _mm_loadu_si128( src + 1 );
			__m128i v2 = _mm_loadu_si128( src + 2 );
			__m128i v3 = _mm_loadu_si128( src + 3 );

			sort4( v0, v1, v2, v3 );

			transpose4x4( v0, v1, v2, v3 );

			_mm_storeu_si128( rdi, v0 );
			_mm_storeu_si128( rdi + 1, v1 );
			_mm_storeu_si128( rdi + 2, v2 );
			_mm_storeu_si128( rdi + 3, v3 );
		}

		// Horizontally sort lanes in the integer vector
		static __forceinline __m128i sortLanes( __m128i v )
		{
			__m128i tmp, i, ax;

			tmp = _mm_shuffle_epi32( v, _MM_SHUFFLE( 1, 0, 3, 2 ) );
			i = Traits::min( v, tmp );
			ax = Traits::max( v, tmp );
			v = _mm_blend_epi16( i, ax, 0b11110000 );

			tmp = _mm_shuffle_epi32( v, _MM_SHUFFLE( 2, 3, 0, 1 ) );
			i = Traits::min( v, tmp );
			ax = Traits::max( v, tmp );
			v = _mm_blend_epi16( i, ax, 0b11001100 );

			tmp = _mm_shuffle_epi32( v, _MM_SHUFFLE( 3, 1, 2, 0 ) );
			i = Traits::min( v, tmp );
			ax = Traits::max( v, tmp );
			return _mm_blend_epi16( i, ax, 0b11110000 );
		}

		// Compare and swap elements of the vector *p1 with the corresponding elements of the vector register *p2, as shown in Figure 5
		static __forceinline void vector_cmpswap( __m128i* p1, __m128i* p2 )
		{
			__m128i a = _mm_loadu_si128( p1 );
			__m128i b = _mm_loadu_si128( p2 );

			__m128i tmp = Traits::min( a, b );
			b = Traits::max( a, b );
			a = tmp;

			_mm_storeu_si128( p1, a );
			_mm_storeu_si128( p2, b );
		}

		// Compare and swap first to third elements of the vector register *p1 with the second to fourth elements of the vector register *p2
		static __forceinline void vector_cmpswap_skew( __m128i* p1, __m128i* p2 )
		{
			__m128i a = _mm_loadu_si128( p1 );
			__m128i b = _mm_loadu_si128( p2 );

			__m128i bs = _mm_srli_si128( b, 4 );
			__m128i a2 = Traits::min( a, bs );
			__m128i b2 = Traits::max( a, bs );
			if constexpr( Traits::minValue == 0 )
			{
				// _mm_srli_si128 shifted in zero, for unsigned integers min( x, 0 ) == x for all x
				a = a2;
			}
			else
			{
				// For negative integers min( x, 0 ) may differ from the x, need to inject the original value
				a = _mm_blend_epi16( a2, a, 0b11000000 );
			}

			b2 = _mm_slli_si128( b2, 4 );
			b = _mm_blend_epi16( b2, b, 0b00000011 );

			_mm_storeu_si128( p1, a );
			_mm_storeu_si128( p2, b );
		}

		// Test whether the array of vectors is vertically sorted
		static __forceinline bool isSorted( const __m128i* rsi, size_t count )
		{
			const __m128i* const rsiEnd = rsi + count;
			__m128i v = _mm_loadu_si128( rsi );
			for( rsi++; rsi < rsiEnd; rsi++ )
			{
				__m128i next = _mm_loadu_si128( rsi );
				__m128i ax = Traits::max( v, next );
				// When values are equal, their XOR is 0
				__m128i xx = _mm_xor_si128( ax, next );
				if( !_mm_testz_si128( xx, xx ) )
					return false;
				v = next;
			}
			return true;
		}

		// Create a vector with broadcasted maximum value for the integers
		static __forceinline __m128i maxVector()
		{
			return _mm_set1_epi32( (int)Traits::maxValue );
		}

		// Load 1-3 elements, pad the remaining ones with maximum value for the scalar
		static __forceinline __m128i loadPartialMaxPad( const void* rsi, size_t length )
		{
			assert( length > 0 && length < 4 );

			__m128i v = maxVector();
			const int* p = (const int*)rsi;
			switch( length )
			{
			case 1:
				return _mm_insert_epi32( v, p[ 0 ], 0 );
			case 2:
				return _mm_insert_epi64( v, *( (int64_t*)p ), 0 );
			case 3:
				v = _mm_insert_epi64( v, *( (int64_t*)p ), 0 );
				return _mm_insert_epi32( v, p[ 2 ], 2 );
			}
			return v;
		}
	};
}  // namespace