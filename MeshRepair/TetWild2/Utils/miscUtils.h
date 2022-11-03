#pragma once

// Copy 3 FP64 values, 24 bytes in total
inline void copyDouble3( double* rdi, const double* rsi )
{
	__m128d v = _mm_loadu_pd( rsi );
	_mm_storeu_pd( rdi, v );

	v = _mm_load_sd( rsi + 2 );
	_mm_store_sd( rdi + 2, v );
}

// Upcast an array of FP32 numbers into FP64. Implemented with SSE2 or AVX1, whichever is enabled in compiler options.
void upcastFloats( double* rdi, size_t length, const float* vb );