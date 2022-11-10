// Load and store routines
#pragma once

namespace AvxMath
{
	// Load 3D vector, set W to 0.0f
	inline __m256d loadDouble3( const double* rsi )
	{
		__m128d low = _mm_loadu_pd( rsi );
		__m128d high = _mm_load_sd( rsi + 2 );
		return _mm256_setr_m128d( low, high );
	}

	// Store 3D vector
	inline void storeDouble3( double* rdi, __m256d vec )
	{
		_mm_storeu_pd( rdi, low2( vec ) );
		_mm_store_sd( rdi + 2, high2( vec ) );
	}
}  // namespace AvxMath