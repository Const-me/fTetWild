#include "stdafx.h"
#include "Box32.h"
using namespace floatTetWild;

namespace
{
	constexpr int fp64MantissaBits = 52;
	constexpr int fp32MantissaBits = 23;

	struct sMiscConstants
	{
		const int signBit = (int)0x80000000u;
		const int maximumPositive = 0x7f800000 - 1;
		const int one = 1;
	};

	alignas( 16 ) static const sMiscConstants miscConstants;

	// Next representable FP32 number, towards +INF
	__forceinline __m128 nextAfter( const __m128 v )
	{
		// FreeBSD source code helped a lot:
		// https://github.com/freebsd/freebsd-src/blob/master/lib/msun/src/s_nextafterf.c

		const __m128i hx = _mm_castps_si128( v );
		const __m128i signBit = _mm_set1_epi32( miscConstants.signBit );
		const __m128i ix = _mm_andnot_si128( signBit, hx );	 // |x|

		__m128i tmp = _mm_set1_epi32( miscConstants.maximumPositive );
		tmp = _mm_cmpgt_epi32( ix, tmp );  // TRUE when the input was INF or NAN

		// Detect INF/NAN, and produce the output for that case
		const __m128 isInfOrNan = _mm_castsi128_ps( tmp );
		const __m128 infOrNanValue = _mm_add_ps( v, v );

		// Detect 0, and produce FLT_MIN to return when it happens
		const __m128 isZero = _mm_cmpeq_ps( v, _mm_setzero_ps() );
		const __m128i oneInt = _mm_set1_epi32( miscConstants.one );
		const __m128 zeroValue = _mm_castsi128_ps( oneInt );

		const __m128 absGreater = _mm_castsi128_ps( _mm_add_epi32( hx, oneInt ) );
		const __m128 absLess = _mm_castsi128_ps( _mm_sub_epi32( hx, oneInt ) );

		__m128 res = _mm_blendv_ps( absGreater, absLess, v );
		res = _mm_blendv_ps( res, zeroValue, isZero );
		res = _mm_blendv_ps( res, infOrNanValue, isInfOrNan );
		return res;
	}

	// Previous representable FP32 number, towards -INF
	__forceinline __m128 prevBefore( const __m128 v )
	{
		// FreeBSD source code helped a lot:
		// https://github.com/freebsd/freebsd-src/blob/master/lib/msun/src/s_nextafterf.c

		const __m128i hx = _mm_castps_si128( v );
		const __m128i signBit = _mm_set1_epi32( miscConstants.signBit );
		const __m128i ix = _mm_andnot_si128( signBit, hx );	 // |x|

		__m128i tmp = _mm_set1_epi32( miscConstants.maximumPositive );
		tmp = _mm_cmpgt_epi32( ix, tmp );  // TRUE when the input was INF or NAN

		// Detect INF/NAN, and produce the output for that case
		const __m128 isInfOrNan = _mm_castsi128_ps( tmp );
		const __m128 infOrNanValue = _mm_add_ps( v, v );

		// Detect 0, and produce -FLT_MIN to return when it happens
		const __m128 isZero = _mm_cmpeq_ps( v, _mm_setzero_ps() );
		const __m128i oneInt = _mm_set1_epi32( miscConstants.one );
		const __m128 zeroValue = _mm_castsi128_ps( _mm_or_si128( signBit, oneInt ) );

		const __m128 absGreater = _mm_castsi128_ps( _mm_add_epi32( hx, oneInt ) );
		const __m128 absLess = _mm_castsi128_ps( _mm_sub_epi32( hx, oneInt ) );

		__m128 res = _mm_blendv_ps( absLess, absGreater, v );
		res = _mm_blendv_ps( res, zeroValue, isZero );
		res = _mm_blendv_ps( res, infOrNanValue, isInfOrNan );
		return res;
	}

	__forceinline __m128 downcastMask( __m256d f64 )
	{
		__m256 f32 = _mm256_castpd_ps( f64 );
		__m128 high = _mm256_extractf128_ps( f32, 1 );
		__m128 low = _mm256_castps256_ps128( f32 );
		return _mm_shuffle_ps( low, high, _MM_SHUFFLE( 2, 0, 2, 0 ) );
	}

	template<bool ceil>
	__forceinline __m128 roundVector( __m256d vec )
	{
		const __m128 rounded = _mm256_cvtpd_ps( vec );
		const __m256d upcasted = _mm256_cvtps_pd( rounded );
		if constexpr( ceil )
		{
			// We want result >= source; compare for upcasted < vec, when true, return the next representable float
			__m256d mask64 = _mm256_cmp_pd( upcasted, vec, _CMP_LT_OQ );
			__m128 next = nextAfter( rounded );
			__m128 mask32 = downcastMask( mask64 );
			return _mm_blendv_ps( rounded, next, mask32 );
		}
		else
		{
			// We want result <= source; compare for upcasted > vec, when true, return the previous representable float
			__m256d mask64 = _mm256_cmp_pd( upcasted, vec, _CMP_GT_OQ );
			__m128 next = prevBefore( rounded );
			__m128 mask32 = downcastMask( mask64 );
			return _mm_blendv_ps( rounded, next, mask32 );
		}
	}

	__forceinline void storeFloat3( float* rdi, __m128 vec )
	{
		_mm_store_sd( (double*)rdi, _mm_castps_pd( vec ) );
		( (int*)rdi )[ 2 ] = _mm_extract_ps( vec, 2 );
	}

	// Make a mask for _mm_insert_ps instruction
	constexpr int insertMask( int source, int dest, int zero = 0 )
	{
		assert( source >= 0 && source < 4 );
		assert( dest >= 0 && dest < 4 );
		assert( zero >= 0 && zero < 16 );
		return ( source << 6 ) | ( dest << 4 ) | zero;
	}

	__forceinline __m128 loadFloat3( const float* rsi )
	{
		__m128 low = _mm_castpd_ps( _mm_load_sd( (const double*)rsi ) );
		__m128 high = _mm_load_ss( rsi + 2 );
		return _mm_insert_ps( low, high, insertMask( 0, 2, 0b1000 ) );
	}
}  // namespace

void Box32::computeTriangle( __m256d a, __m256d b, __m256d c )
{
	const __m256d i = _mm256_min_pd( _mm256_min_pd( a, b ), c );
	const __m256d ax = _mm256_max_pd( _mm256_max_pd( a, b ), c );

	__m128 i32 = roundVector<false>( i );
	__m128 ax32 = roundVector<true>( ax );
	_mm_storeu_ps( &boxMin[ 0 ], i32 );
	storeFloat3( &boxMax[ 0 ], ax32 );
}

void Box32::computeUnion( const Box32& a, const Box32& b )
{
	__m128 v1 = _mm_loadu_ps( &a.boxMin[ 0 ] );
	__m128 v2 = _mm_loadu_ps( &b.boxMin[ 0 ] );
	_mm_storeu_ps( &boxMin[ 0 ], _mm_min_ps( v1, v2 ) );

	v1 = loadFloat3( &a.boxMax[ 0 ] );
	v2 = loadFloat3( &b.boxMax[ 0 ] );
	storeFloat3( &boxMax[ 0 ], _mm_max_ps( v1, v2 ) );
}

namespace
{
	// [ x, y, z, w ] => [ min( x, y, z ), y, z, w ]
	__forceinline __m128 horizontalMinimum3( __m128 v )
	{
		__m128 z = _mm_movehl_ps( v, v );
		__m128 y = _mm_movehdup_ps( v );
		v = _mm_min_ss( v, z );
		v = _mm_min_ss( v, y );
		return v;
	}

	// Compute bitwise OR of all 4 lanes in the vector, broadcast the result
	__forceinline __m128 horizontalBitwiseOr( __m128 v )
	{
		// v |= v.yxwz
		__m128 p = _mm_permute_ps( v, _MM_SHUFFLE( 2, 3, 0, 1 ) );
		v = _mm_or_ps( v, p );

		// v |= v.zwxy
		p = _mm_permute_ps( v, _MM_SHUFFLE( 1, 0, 3, 2 ) );
		v = _mm_or_ps( v, p );
		return v;
	}
}  // namespace

__m128 Box32::pointBoxSignedSquaredDistance( __m128 pos ) const
{
	// Load the box into 2 vectors
	const __m128 boxMin = _mm_loadu_ps( &this->boxMin[ 0 ] );
	const __m128 boxMax = _mm_loadu_ps( &this->boxMax[ 0 ] );

	// When inside, both numbers are positive
	const __m128 dmin = _mm_sub_ps( pos, boxMin );
	const __m128 dmax = _mm_sub_ps( boxMax, pos );

	// When inside, distance to the box
	// When outside, one of the vectors was negative another one positive, min will return the negative one, which is the distance we're after
	__m128 dist = _mm_min_ps( dmin, dmax );

	// Compute mask of the lanes which were outside the box
	const __m128 zero = _mm_setzero_ps();
	const __m128 outsideMaskVector = _mm_cmplt_ps( dist, zero );

	// ---- Compute positive distance for the outside case ----
	__m128 resOut = _mm_and_ps( dist, outsideMaskVector );
	resOut = _mm_dp_ps( resOut, resOut, 0b01111111 );

	// ---- Compute negative distance for the inside case ----
	__m128 resIn = horizontalMinimum3( dist );
	// Compute square
	resIn = _mm_mul_ss( resIn, resIn );
	// Negate the vector
	resIn = _mm_sub_ss( zero, resIn );
	// Broadcast the lower lane across the complete vector
#ifdef __AVX2__
	resIn = _mm_broadcastss_ps( resIn );
#else
	resIn = _mm_permute_ps( resIn, _MM_SHUFFLE( 0, 0, 0, 0 ) );
#endif

	// ---- Make a mask to select between the two ----
	__m128 mask = _mm_blend_ps( outsideMaskVector, zero, 0b1000 );
	mask = horizontalBitwiseOr( mask );

	// Select one of the output vectors
	return _mm_blendv_ps( resIn, resOut, mask );
}

namespace
{
	__forceinline __m256 horizontalMinimum3( __m256 v )
	{
		__m256 z = _mm256_permute_ps( v, _MM_SHUFFLE( 3, 2, 1, 2 ) );
		__m256 y = _mm256_movehdup_ps( v );
		v = _mm256_min_ps( v, z );
		v = _mm256_min_ps( v, y );
		return v;
	}

	// Compute bitwise OR of 4 lanes in the vector in each half, broadcast the result across 16-byte lanes
	__forceinline __m256 horizontalBitwiseOr( __m256 v )
	{
		// v |= v.yxwz
		__m256 p = _mm256_permute_ps( v, _MM_SHUFFLE( 2, 3, 0, 1 ) );
		v = _mm256_or_ps( v, p );

		// v |= v.zwxy
		p = _mm256_permute_ps( v, _MM_SHUFFLE( 1, 0, 3, 2 ) );
		v = _mm256_or_ps( v, p );
		return v;
	}

	// Slightly faster equivalent of _mm256_dp_ps( v, v, 0b01111111 )
	__forceinline __m256 len2( __m256 v )
	{
#if 1
		// On Zen 3, vdpps has 15 cycles latency and 7 micro-ops
		// vmulps = 3 cycles, vpermilps = 3 cycles, vaddps = 3 cycles, so this version is hopefully 12 cycles combined, and 6 micro-ops
		v = _mm256_mul_ps( v, v );
		const __m256 x = _mm256_permute_ps( v, _MM_SHUFFLE( 0, 0, 0, 0 ) );
		const __m256 y = _mm256_permute_ps( v, _MM_SHUFFLE( 1, 1, 1, 1 ) );
		const __m256 z = _mm256_permute_ps( v, _MM_SHUFFLE( 2, 2, 2, 2 ) );
		const __m256 xy = _mm256_add_ps( x, y );
		return _mm256_add_ps( xy, z );
#else
		// If your computer has Intel CPU instead of AMD, try to use this version instead, and profile.
		return _mm256_dp_ps( v, v, 0b01111111 );
#endif
	}
}  // namespace

__m256 Box32::pointBoxSignedSquaredDistanceX2( const Box32& b1, const Box32& b2, __m256 pos )
{
	// Load both boxes into 2 AVX vectors
	const __m128 boxMin1 = _mm_loadu_ps( &b1.boxMin[ 0 ] );
	const __m128 boxMax1 = _mm_loadu_ps( &b1.boxMax[ 0 ] );
	const __m128 boxMin2 = _mm_loadu_ps( &b2.boxMin[ 0 ] );
	const __m128 boxMax2 = _mm_loadu_ps( &b2.boxMax[ 0 ] );
	const __m256 boxMin = _mm256_setr_m128( boxMin1, boxMin2 );
	const __m256 boxMax = _mm256_setr_m128( boxMax1, boxMax2 );

	// When inside, both numbers are positive
	const __m256 dmin = _mm256_sub_ps( pos, boxMin );
	const __m256 dmax = _mm256_sub_ps( boxMax, pos );

	// When inside, distance to the box
	// When outside, one of the vectors was negative another one positive, min will return the negative one, which is the distance we're after
	__m256 dist = _mm256_min_ps( dmin, dmax );

	// Compute mask of the lanes which were outside the box
	const __m256 zero = _mm256_setzero_ps();
	const __m256 outsideMaskVector = _mm256_cmp_ps( dist, zero, _CMP_LT_OQ );

	// ---- Compute positive distance for the outside case ----
	__m256 resOut = _mm256_and_ps( dist, outsideMaskVector );
	resOut = len2( resOut );

	// ---- Compute negative distance for the inside case ----
	__m256 resIn = horizontalMinimum3( dist );
	// Compute square
	resIn = _mm256_mul_ps( resIn, resIn );
	// Negate the vector
	resIn = _mm256_sub_ps( zero, resIn );
	// Broadcast the lower lane across the complete vector
	resIn = _mm256_permute_ps( resIn, _MM_SHUFFLE( 0, 0, 0, 0 ) );

	// ---- Make a mask to select between the two ----
	__m256 mask = _mm256_blend_ps( outsideMaskVector, zero, 0b10001000 );
	mask = horizontalBitwiseOr( mask );

	// Select one of the output vectors
	return _mm256_blendv_ps( resIn, resOut, mask );
}