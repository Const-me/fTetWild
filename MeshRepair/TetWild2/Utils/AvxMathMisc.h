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

	// Compute horizontal sum of the 3D vector
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

	// Compute square of the distance between two 3D vectors
	inline double vector3DistanceSquared( __m256d a, __m256d b )
	{
		const __m256d diff = _mm256_sub_pd( b, a );
		return vector3DotScalar( diff, diff );
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

	// Compute two pieces of the cross-product; the complete cross product is equal to ( lhs - rhs )
	__forceinline void crossProductPieces( __m256d a, __m256d b, __m256d& lhs, __m256d& rhs )
	{
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
		lhs = _mm256_mul_pd( a1, b2 );
		rhs = _mm256_mul_pd( a2, b1 );
	}

	// Compute cross product between two 3D vectors
	// The unused W lane is set to 0 unless there's INF or NAN in W lanes of the inputs
	inline __m256d vector3Cross( __m256d a, __m256d b )
	{
		// The formula is a.yzx * b.zxy - a.zxy * b.yzx
		// https://en.wikipedia.org/wiki/Cross_product#Coordinate_notation
		__m256d lhs, rhs;
		crossProductPieces( a, b, lhs, rhs );
		return _mm256_sub_pd( lhs, rhs );
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

	// Performs a linear interpolation using the incorrect but fast formula:  x * ( 1 - s ) + y * s
	inline __m256d lerpFast( __m256d x, __m256d y, __m256d s )
	{
		__m256d s0 = _mm256_sub_pd( _mm256_set1_pd( 1.0 ), s );
		x = _mm256_mul_pd( x, s0 );
		y = _mm256_mul_pd( y, s );
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
#ifdef __AVX2__

		// _mm256_permutevar_pd only shuffles within 16-byte pieces, need _mm256_permutevar8x32_ps
		// Create int32 permutation vector, first 2 lanes equal to [ lane * 2, lane * 2 + 1 ]
		uint64_t perm = lane * 2;
		perm |= ( ( perm | 1 ) << 32 );
		__m128i v1 = _mm_cvtsi64_si128( (int64_t)perm );
		// We don't care what's in the other 6 lanes, undefined is fine
		__m256i v2 = _mm256_castsi128_si256( v1 );

		// Cast to FP32 vector, and shuffle with that permutation
		__m256 f = _mm256_castpd_ps( vec );
		f = _mm256_permutevar8x32_ps( f, v2 );
		// Cast back to FP64 vector, and return the lowest lane
		vec = _mm256_castps_pd( f );
		return _mm256_cvtsd_f64( vec );
#else
		// Select either high or low slice of the input vector
		const __m128d high = _mm256_extractf128_pd( vec, 1 );
		const __m128d low = _mm256_castpd256_pd128( vec );

		const bool selectBool = ( 0 != ( lane & 2 ) );
		const __m128i selectMask = _mm_set1_epi32( selectBool ? -1 : 0 );
		const __m128d vec2 = _mm_blendv_pd( low, high, _mm_castsi128_pd( selectMask ) );

		// Use _mm_permutevar_pd instruction to return either x or y lane of the vec2
		// Funfact: just for this one instruction Intel ignores the lowest bit of the second argument, and the index starts at bit #1
		// Probably a hardware bug in Sandy Bridge (first AVX1 CPU) discovered too late, so Intel documented the bug and called it a day.
		const __m128i v1 = _mm_cvtsi32_si128( ( lane & 1 ) << 1 );
		const __m128d res = _mm_permutevar_pd( vec2, v1 );
		return _mm_cvtsd_f64( res );
#endif
	}

	inline __m256d vector3Normalize( __m256d v )
	{
		__m256d div = _mm256_set1_pd( vector3Length( v ) );
		return _mm256_div_pd( v, div );
	}

	// Compute ( a * b ) + c, using FMA in AVX2 builds
	inline __m256d vectorMultiplyAdd( __m256d a, __m256d b, __m256d c )
	{
#ifdef __AVX2__
		return _mm256_fmadd_pd( a, b, c );
#else
		return _mm256_add_pd( _mm256_mul_pd( a, b ), c );
#endif
	}

	// Extract X lane from the vector
	inline double vectorGetX( __m256d vec )
	{
		return _mm256_cvtsd_f64( vec );
	}

	// Extract Y lane from the vector
	inline double vectorGetY( __m256d vec )
	{
		return _mm_cvtsd_f64( _mm_permute_pd( low2( vec ), 0b11 ) );
	}

	// Extract Z lane from the vector
	inline double vectorGetZ( __m256d vec )
	{
		return _mm_cvtsd_f64( high2( vec ) );
	}

	// Extract W lane from the vector
	inline double vectorGetW( __m256d vec )
	{
		__m128d high = high2( vec );
		high = _mm_permute_pd( high, 0b11 );
		return _mm_cvtsd_f64( high );
	}

	// Compute squared lengths of 3 vectors
	inline __m256d computeLengthsSquared( __m256d v0, __m256d v1, __m256d v2 )
	{
		// Reduce count of shuffles by transposing them first
		__m256d c0, c1, c2;
		transpose3x3( v0, v1, v2, c0, c1, c2 );
		// Once transposed, we need vertical operations instead of horizontal, faster
		c0 = _mm256_mul_pd( c0, c0 );
		c1 = _mm256_mul_pd( c1, c1 );
		c2 = _mm256_mul_pd( c2, c2 );

		return _mm256_add_pd( _mm256_add_pd( c0, c1 ), c2 );
	}

	// Compute lengths of 3 vectors
	inline __m256d computeLengths( __m256d v0, __m256d v1, __m256d v2 )
	{
		const __m256d len2 = computeLengthsSquared( v0, v1, v2 );
		return _mm256_sqrt_pd( len2 );
	}
}  // namespace AvxMath