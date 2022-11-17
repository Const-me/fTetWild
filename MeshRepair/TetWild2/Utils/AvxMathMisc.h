#pragma once
#include <assert.h>
#include <immintrin.h>

namespace AvxMath
{
	// Extract [ X, Y ] slice of the vector; the function is free in runtime, compiles into no instructions
	inline __m128d low2( __m256d vec )
	{
		return _mm256_castpd256_pd128( vec );
	}

	// Extract [ Z, W ] slice of the vector
	inline __m128d high2( __m256d vec )
	{
		return _mm256_extractf128_pd( vec, 1 );
	}

	inline double vector3HorizontalSum( __m256d v )
	{
		__m128d xy = low2( v );
		__m128d z = high2( v );
		xy = _mm_add_sd( xy, _mm_unpackhi_pd( xy, xy ) );
		xy = _mm_add_sd( xy, z );
		return _mm_cvtsd_f64( xy );
	}

	// Horizontal sum of the 3D vector, broadcasted to both lanes of SSE vector
	inline __m128d vector3HorizontalSum2( __m256d v )
	{
		__m128d xy = low2( v );
		__m128d z = high2( v );
		xy = _mm_add_sd( xy, _mm_unpackhi_pd( xy, xy ) );
		xy = _mm_add_sd( xy, z );
		return _mm_movedup_pd( xy );
	}

	// Horizontal sum of the 3D vector, return an SSE vector
	inline __m128d vector3HorizontalSum1( __m256d v )
	{
		__m128d xy = low2( v );
		__m128d z = high2( v );
		xy = _mm_add_sd( xy, _mm_unpackhi_pd( xy, xy ) );
		xy = _mm_add_sd( xy, z );
		return xy;
	}

	// Compute dot product of two 3D vectors, return a scalar
	inline double vector3DotScalar( __m256d a, __m256d b )
	{
		const __m256d prod = _mm256_mul_pd( a, b );
		return vector3HorizontalSum( prod );
	}

	inline double vector3Length( __m256d v )
	{
		__m256d prod = _mm256_mul_pd( v, v );
		__m128d hs = vector3HorizontalSum1( prod );
		__m128d res = _mm_sqrt_sd( hs, hs );
		return _mm_cvtsd_f64( res );
	}

	// Dot product of 3D vectors, broadcast to both lanes of SSE vector
	inline __m128d vector3Dot2( __m256d a, __m256d b )
	{
		const __m256d prod = _mm256_mul_pd( a, b );
		return vector3HorizontalSum2( prod );
	}

	// Permute XYZW => ZWXY
	inline __m256d flipHighLow( __m256d vec )
	{
		return _mm256_permute2f128_pd( vec, vec, 1 );
	}

	// Transpose 3x3 FP64 matrix in 3 registers. The 4-th W lane of the result is garbage.
	__forceinline void transpose3x3( __m256d a, __m256d b, __m256d c, __m256d& x, __m256d& y, __m256d& z )
	{
		__m256d axbxazbz = _mm256_unpacklo_pd( a, b );			 // a.x, b.x, a.z, b.z
		z = _mm256_permute2f128_pd( axbxazbz, c, 0x31 );		 // a.z, b.z, c.z, c.w
		__m128d ayby = _mm_unpackhi_pd( low2( a ), low2( b ) );	 // a.y, b.y
		x = _mm256_insertf128_pd( axbxazbz, low2( c ), 1 );		 // a.x, b.x, c.x, c.y
		__m128d cyy = _mm_unpackhi_pd( low2( c ), low2( c ) );	 // c.y, c.y
		y = _mm256_setr_m128d( ayby, cyy );						 // a.y, b.y, cy, cy
	}

	// Compute cross product between two 3D vectors
	// The unused W lane is set to 0 unless there's INF or NAN in W lanes of the inputs
	inline __m256d vector3Cross( __m256d a, __m256d b )
	{
		// The formula is a.yzx * b.zxy - a.zxy * b.yzx
		// https://en.wikipedia.org/wiki/Cross_product#Coordinate_notation

#ifdef __AVX2__
		const __m256d a1 = _mm256_permute4x64_pd( a, _MM_SHUFFLE( 3, 0, 2, 1 ) );  // a.yzxw
		const __m256d b2 = _mm256_permute4x64_pd( b, _MM_SHUFFLE( 3, 1, 0, 2 ) );  // b.zxyw
		const __m256d a2 = _mm256_permute4x64_pd( a, _MM_SHUFFLE( 3, 1, 0, 2 ) );  // a.zxyw
		const __m256d b1 = _mm256_permute4x64_pd( b, _MM_SHUFFLE( 3, 0, 2, 1 ) );  // b.yzxw
#else
		const __m256d af = flipHighLow( a );  // a.zwxy
		const __m256d bf = flipHighLow( b );  // b.zwxy

		__m256d a1 = _mm256_shuffle_pd( a, af, 0b0101 );		// a.yzwx
		__m256d b1 = _mm256_shuffle_pd( b, bf, 0b0101 );		// b.yzwx
		const __m256d b2 = _mm256_shuffle_pd( bf, b, 0b1100 );	// b.zxyw
		const __m256d a2 = _mm256_shuffle_pd( af, a, 0b1100 );	// a.zxyw
		a1 = _mm256_permute_pd( a1, 0b0110 );					// a.yzxw
		b1 = _mm256_permute_pd( b1, 0b0110 );					// b.yzxw
#endif
		return _mm256_sub_pd( _mm256_mul_pd( a1, b2 ), _mm256_mul_pd( a2, b1 ) );
	}

	// Negate all 4 lanes of the vector
	inline __m256d vectorNegate( __m256d v )
	{
		__m256d zero = _mm256_setzero_pd();
		return _mm256_sub_pd( zero, v );
	}

	// Absolute value of all 4 lanes
	inline __m256d vectorAbs( __m256d v )
	{
		__m256d neg = vectorNegate( v );
		return _mm256_max_pd( v, neg );
	}

	// Performs a linear interpolation using the incorrect but fast formula:  x * ( 1 - s ) + y * s
	// The correct formula is written in that document: https://open-std.org/jtc1/sc22/wg21/docs/papers/2019/p0811r3.html
	inline __m256d lerpFast( __m256d x, __m256d y, double s )
	{
		x = _mm256_mul_pd( x, _mm256_set1_pd( 1.0 - s ) );
		y = _mm256_mul_pd( y, _mm256_set1_pd( s ) );
		return _mm256_add_pd( x, y );
	}

	inline __m256d vectorBroadcast( __m128d x )
	{
#ifdef __AVX2__
		return _mm256_broadcastsd_pd( x );
#else
		x = _mm_movedup_pd( x );
		return _mm256_setr_m128d( x, x );
#endif
	}

	// Compute maximum of the 3 XYZ lanes, and broadcast the value over all 4 lanes of the result
	inline __m256d vector3BroadcastMaximum( __m256d v )
	{
		// Compute maximum of XYZ lanes
		const __m128d high = _mm256_extractf128_pd( v, 1 );
		const __m128d low = _mm256_castpd256_pd128( v );
		__m128d max2 = _mm_max_sd( low, _mm_unpackhi_pd( low, low ) );
		max2 = _mm_max_sd( max2, high );
		// Broadcast the maximum to all lanes of another vector
		return vectorBroadcast( max2 );
	}

	// Compare vectors for a == b, return index of the first lane which compared as true.
	// If none of the lanes were equal, returns 32.
	inline uint32_t firstEqualLaneIndex( __m256d a, __m256d b )
	{
		// Compare for equality with n
		const __m256d eq = _mm256_cmp_pd( a, b, _CMP_EQ_OQ );
		// Return index of the first equal lane
		const uint32_t mask = _mm256_movemask_pd( eq );
		return (int)_tzcnt_u32( mask );
	}

	// Extract specified lane from the vector
	inline double vectorExtractLane( __m256d vec, uint32_t lane )
	{
		assert( lane < 4 );
		__m128i v1 = _mm_cvtsi32_si128( (int)lane );
		__m256i v2 = _mm256_castsi128_si256( v1 );
		// Surprisingly, that instruction is from AVX1 set, no need for #ifdef-s
		vec = _mm256_permutevar_pd( vec, v2 );
		return _mm256_cvtsd_f64( vec );
	}

	inline __m256d vector3Normalize( __m256d v )
	{
		__m256d div = _mm256_set1_pd( vector3Length( v ) );
		return _mm256_div_pd( v, div );
	}
}  // namespace AvxMath