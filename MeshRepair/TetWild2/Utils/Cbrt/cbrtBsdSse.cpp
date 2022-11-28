#include "stdafx.h"
#include "cbrt.h"
#include <stdint.h>
#include <immintrin.h>

namespace
{
	constexpr uint32_t B1 = 715094163;	// B1 = (1023-1023/3-0.03306235651)*2**20
	constexpr uint32_t B2 = 696219795;	// B2 = (1023-1023/3-54/3-0.03306235651)*2**20

	// Get the more significant 32 bit int from the lower FP64 lane of the vector
	__forceinline uint32_t GET_HIGH_WORD( __m128d vec )
	{
		__m128i iv = _mm_castpd_si128( vec );
		return (uint32_t)_mm_extract_epi32( iv, 1 );
	}

	// Some magic numbers used in the process
	// Reordering these constants is guaranteed to break the implementation. Most of them are loaded two at a time, with 128-bit vector loads
	struct sCbrtConstants
	{
		const double P3 = -0.758397934778766047437;			   // 0xbfe844cb, 0xbee751d9
		const double P1 = -1.88497979543377169875;			   // 0xbffe28e0, 0x92f02420
		const double P4 = 0.145996192886612446982;			   // 0x3fc2b000, 0xd4e4edd7
		const double P2 = 1.621429720105354466140;			   // 0x3ff9f160, 0x4a49d6c2
		const double P0 = 1.87595182427177009643;			   // 0x3ffe03e6, 0x0f61e692
		const uint64_t denormMagic = ( 0x43500000ull << 32 );  // 2**54
	};

	// Without two extra vectors all these constants fit in a single cache line, with 16 unused bytes.
	alignas( 64 ) static const sCbrtConstants cbrtConstants;
}  // namespace


// ====================================================
// Copyright (c) 1993 by Sun Microsystems, Inc. All rights reserved.
// Developed at SunPro, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this software is freely granted, provided that this notice is preserved.
// Optimized by Bruce D.Evans; then optimized by Konstantin, http://const.me
double __declspec( noinline ) LowLevel::cbrtBsdSse( double val )
{
	const __m128d x = _mm_set_sd( val );
	__m128i vi = _mm_castpd_si128( x );

	uint32_t hx = _mm_extract_epi32( vi, 1 );
	constexpr uint32_t signBit = 0x80000000u;
	const uint32_t sign = hx & signBit;

	// Assuming BMI1 support: https://en.wikipedia.org/wiki/X86_Bit_manipulation_instruction_set#Supporting_CPUs
	hx = _andn_u32( signBit, hx );

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
	uint32_t offsetNumber;
	vi = _mm_setzero_si128();
	if( hx >= 0x00100000 )
		offsetNumber = B1;
	else
	{
		if( _mm_comieq_sd( x, _mm_setzero_pd() ) )
			return _mm_cvtsd_f64( x );	// cbrt(0) is itself
		// set t = 2**54
		t = _mm_load_sd( (const double*)&cbrtConstants.denormMagic );
		t = _mm_mul_sd( x, t );
		hx = GET_HIGH_WORD( t );
		hx = _andn_u32( signBit, hx );
		offsetNumber = B2;
	}
	vi = _mm_insert_epi32( vi, (int)( sign | ( hx / 3 + offsetNumber ) ), 1 );
	t = _mm_castsi128_pd( vi );

	// New cbrt to 23 bits:
	//     cbrt(x) = t*cbrt(x/t**3) ~= t*P(t**3/x)
	// where P(r) is a polynomial of degree 4 that approximates 1/cbrt(r) to within 2**-23.5 when |r - 1| < 1/10.
	// The rough approximation has produced t such than |t/cbrt(x) - 1| ~< 1/32, and cubing this gives us bounds for r = t**3/x.

	// r = ( t * t ) * ( t / x );
	__m128d r = _mm_mul_sd( _mm_mul_sd( t, t ), _mm_div_sd( t, x ) );

	// Original code:
	// t = t * ( ( P0 + r * ( P1 + r * P2 ) ) + ( ( r * r ) * r ) * ( P3 + r * P4 ) );
	const __m128d rr = _mm_movedup_pd( r );				   // [ r, r ]
	const __m128d r2r = _mm_mul_sd( rr, rr );			   // [ r*r, r ]
	const __m128d r3r = _mm_mul_sd( rr, r2r );			   // [ r*r*r, r ]
	const __m128d p31 = _mm_load_pd( &cbrtConstants.P3 );  // [ P3, P1 ]
	const __m128d p42 = _mm_load_pd( &cbrtConstants.P4 );  // [ P4, P2 ]
	// We don't really want FMA
	// We want compatibility with the original code from BSD
	// Also, with two separate add & multiply instructions both memory loads are merged, with FMA we would need a separate load instruction
	__m128d rhs = _mm_add_pd( _mm_mul_pd( rr, p42 ), p31 );	 // [ P3 + r * P4, P1 + r * P2 ]

	rhs = _mm_mul_pd( r3r, rhs );								// [ r*r*r * ( P3 + r * P4), r*( P1 + r * P2 ) ]
	__m128d tmp = _mm_unpackhi_pd( rhs, rhs );					// r*( P1 + r * P2 )
	tmp = _mm_add_sd( tmp, _mm_load_sd( &cbrtConstants.P0 ) );	// P0 + higher lane of the rhs
	tmp = _mm_add_sd( tmp, rhs );								// Add to the low lane of the rhs
	t = _mm_mul_sd( t, tmp );

	// Round t away from zero to 23 bits (sloppily except for ensuring that the result is larger in magnitude than cbrt(x) but not much more than 2 23-bit
	// ulps larger). With rounding towards zero, the error bound would be ~5/6 instead of ~4/6. With a maximum error of 2 23-bit ulps in the rounded t, the
	// infinite-precision error in the Newton approximation barely affects third digit in the final error 0.667; the error in the rounded t can be up to
	// about 3 23-bit ulps before the final error is larger than 0.667 ulps.

	// Original code: u.bits = ( u.bits + 0x80000000 ) & 0xffffffffc0000000ULL;
	vi = _mm_castpd_si128( t );
	uint64_t bits = (uint64_t)_mm_cvtsi128_si64( vi );
	bits += 0x80000000;	// Same value as signBit couple pages above, VC++ is smart enough to keep that number in a register
	bits &= 0xFFFFFFFFC0000000ull;
	vi = _mm_cvtsi64_si128( (int64_t)bits );
	t = _mm_castsi128_pd( vi );

	// One step Newton iteration to 53 bits with error < 0.667 ulps
	// The remaining code is scalar, not sure whether it's possible to optimize somehow
	double ts = _mm_cvtsd_f64( t );
	double xs = _mm_cvtsd_f64( x );

	double s = ts * ts;				// t * t is exact
	double rs = xs / s;				// error <= 0.5 ulps; |r| < |t|
	double w = ts + ts;				// t+t is exact
	rs = ( rs - ts ) / ( w + rs );	// r-t is exact; w+r ~= 3*t
	ts = ts + ts * rs;				// error <= 0.5 + 0.5/3 + epsilon
	return ts;
}