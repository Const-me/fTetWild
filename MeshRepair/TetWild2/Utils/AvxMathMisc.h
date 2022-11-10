#pragma once

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

	// Compute dot product of two 3D vectors, return a scalar
	inline double vector3DotScalar( __m256d a, __m256d b )
	{
		const __m256d prod = _mm256_mul_pd( a, b );
		__m128d xy = low2( prod );
		__m128d z = high2( prod );
		xy = _mm_add_sd( xy, _mm_unpackhi_pd( xy, xy ) );
		xy = _mm_add_sd( xy, z );
		return _mm_cvtsd_f64( xy );
	}
}  // namespace AvxMath