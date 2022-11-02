#pragma once

inline void copyDouble3( double* rdi, const double* rsi )
{
	__m128d v = _mm_loadu_pd( rsi );
	_mm_storeu_pd( rdi, v );

	v = _mm_load_sd( rsi + 2 );
	_mm_store_sd( rdi + 2, v );
}