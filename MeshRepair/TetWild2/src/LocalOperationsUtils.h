#pragma once

namespace floatTetWild
{
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
}