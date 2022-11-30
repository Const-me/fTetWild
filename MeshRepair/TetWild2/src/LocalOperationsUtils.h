#pragma once
#include <Utils/AvxMathVec.h>

namespace floatTetWild
{
	using namespace AvxMath;

	__forceinline void translateToBarycenter( __m256d& v0, __m256d& v1, __m256d& v2, __m256d& v3, __m256d oneFourth )
	{
		// Compute barycenter of the vertices
		const __m256d s01 = _mm256_add_pd( v0, v1 );
		const __m256d s23 = _mm256_add_pd( v2, v3 );
		const __m256d s = _mm256_add_pd( s01, s23 );
		const __m256d bc = _mm256_mul_pd( s, oneFourth );

		// Translate the input vertices, making them relative to the barycenter
		v0 = _mm256_sub_pd( v0, bc );
		v1 = _mm256_sub_pd( v1, bc );
		v2 = _mm256_sub_pd( v2, bc );
		v3 = _mm256_sub_pd( v3, bc );
	}

	// Broadcast a scalar from memory to all lanes of the vector
	inline __m256d broadcast( const double& r )
	{
		return _mm256_broadcast_sd( &r );
	}

	// Add 12 numbers in XYZ lanes of 4 vectors
	inline double hadd12( __m256d a, __m256d b, __m256d c, __m256d d )
	{
		// Compute vertical sum first
		const __m256d ab = _mm256_add_pd( a, b );
		const __m256d cd = _mm256_add_pd( c, d );
		const __m256d v = _mm256_add_pd( ab, cd );
		// Now add the 3 numbers
		return vector3HorizontalSum( v );
	}

	// Compute x - y + z
	inline double horizontalAddSub( __m256d v )
	{
		__m128d z = _mm256_extractf128_pd( v, 1 );
		__m128d xy = _mm256_castpd256_pd128( v );
		xy = _mm_sub_sd( xy, _mm_unpackhi_pd( xy, xy ) );
		xy = _mm_add_sd( xy, z );
		return _mm_cvtsd_f64( xy );
	}

#ifdef __AVX2__
	inline __m256d permute_yxx( __m256d vec )
	{
		return _mm256_permute4x64_pd( vec, _MM_SHUFFLE( 3, 0, 0, 1 ) );	 // yxxw
	}
	inline __m256d permute_zzy( __m256d vec )
	{
		return _mm256_permute4x64_pd( vec, _MM_SHUFFLE( 3, 1, 2, 2 ) );	 // zzyw
	}
#else
	inline __m256d permute_yxx( __m256d vec )
	{
		const __m128d xy = _mm256_castpd256_pd128( vec );
		const __m128d yx = _mm_permute_pd( xy, _MM_SHUFFLE2( 0, 1 ) );
		return _mm256_setr_m128d( yx, xy );	 // yxxy
	}
	inline __m256d permute_zzy( __m256d vec )
	{
		const __m256d zwxy = _mm256_permute2f128_pd( vec, vec, 1 );
		return _mm256_permute_pd( zwxy, 0b1100 );  // zzyy
	}
#endif	// __AVX2__

	inline __m256d customProduct( __m256d a, __m256d b )
	{
		// The normal cross product formula is a.yzx * b.zxy - a.zxy * b.yzx
		// These mysterios AMIPS formulae instead computing a.yxx * b.zzy - a.zzy * b.yxx, systematically so in quite a few places
		// This results in cross product with inverted Y axis
		const __m256d a0 = permute_yxx( a );
		const __m256d b0 = permute_zzy( b );
		const __m256d a1 = permute_zzy( a );
		const __m256d b1 = permute_yxx( b );
		const __m256d cp1 = _mm256_mul_pd( a0, b0 );
		const __m256d cp2 = _mm256_mul_pd( a1, b1 );
		return _mm256_sub_pd( cp1, cp2 );
	}

	// Extract XYZ lanes from the vector, and store into 3 separate const variables with _x, _y and _z suffix
#define STORE( v )                        \
	const double v##_x = vectorGetX( v ); \
	const double v##_y = vectorGetY( v ); \
	const double v##_z = vectorGetZ( v );

}  // namespace floatTetWild