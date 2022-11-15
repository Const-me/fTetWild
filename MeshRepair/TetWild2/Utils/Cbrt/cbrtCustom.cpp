#include "stdafx.h"
#include "cbrt.h"
#include <stdint.h>
#include <immintrin.h>

namespace
{
	constexpr uint32_t B1 = 715094163;	// B1 = (1023-1023/3-0.03306235651)*2**20
	constexpr uint32_t B2 = 696219795;	// B2 = (1023-1023/3-54/3-0.03306235651)*2**20

	// Get the more significant 32 bit int from a double.
	__forceinline uint32_t GET_HIGH_WORD( __m128d vec )
	{
		__m128i iv = _mm_castpd_si128( vec );
		return (uint32_t)_mm_extract_epi32( iv, 1 );
	}

	__forceinline __m128d setIntegerFp64( uint64_t val )
	{
		__m128i v = _mm_cvtsi64_si128( (int64_t)val );
		return _mm_castsi128_pd( v );
	}

	// iterative cube root approximation using Halley's method (double)
	__forceinline double cbrta_halleyd( const double a, const double R )
	{
		const double a3 = a * a * a;
		const double b = a * ( a3 + R + R ) / ( a3 + a3 + R );
		return b;
	}
}  // namespace

double LowLevel::cbrtCustom( double val )
{
	__m128d x = _mm_set_sd( val );
	__m128i vi = _mm_castpd_si128( x );

	int hx = _mm_extract_epi32( vi, 1 );
	constexpr uint32_t signBit = 0x80000000u;
	const uint32_t sign = hx & signBit;

#ifdef __AVX__
	// Assuming BMI1 support: https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Supporting_CPUs
	hx = (int)_andn_u32( signBit, (uint32_t)hx );
#else
	hx &= ( ~signBit );
#endif

	if( hx >= 0x7ff00000 )
	{
		// cbrt(NaN,INF) is itself
		return _mm_cvtsd_f64( _mm_add_sd( x, x ) );
	}

	__m128d t;
	// Rough cbrt to 5 bits:
	// cbrt(2**e*(1+m) ~= 2**(e/3)*(1+(e%3+m)/3)
	// where e is integral and >= 0, m is real and in [0, 1), and "/" and "%" are integer division and modulus with rounding towards minus infinity.
	// The RHS is always >= the LHS and has a maximum relative error of about 1 in 16.
	// Adding a bias of -0.03306235651 to the (e%3+m)/3 term reduces the error to about 1 in 32.
	// With the IEEE floating point representation, for finite positive normal values, ordinary integer division of the value in bits
	// magically gives almost exactly the RHS of the above provided we first subtract the exponent bias (1023 for doubles) and later add it back.
	// We do the subtraction virtually to keep e >= 0 so that ordinary integer division rounds towards minus infinity; this is also efficient.
	vi = _mm_setzero_si128();
	if( hx < 0x00100000 )
	{
		if( _mm_comieq_sd( x, _mm_setzero_pd() ) )
			return _mm_cvtsd_f64( x );	// cbrt(0) is itself

		// set t = 2**54
		t = setIntegerFp64( 0x4350000000000000ull );
		t = _mm_mul_sd( x, t );
		uint32_t high = GET_HIGH_WORD( t );
		high = _andn_u32( signBit, high );
		vi = _mm_insert_epi32( vi, (int)( sign | ( high / 3 + B2 ) ), 1 );
	}
	else
	{
		vi = _mm_insert_epi32( vi, (int)( sign | ( hx / 3 + B1 ) ), 1 );
	}

	t = _mm_castsi128_pd( vi );
	const double R = _mm_cvtsd_f64( x );
	double a = _mm_cvtsd_f64( t );

	a = cbrta_halleyd( a, R );
	a = cbrta_halleyd( a, R );
	return cbrta_halleyd( a, R );
}