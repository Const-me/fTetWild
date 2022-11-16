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

#ifdef __AVX__
// Compute horizontal sum of FP64 elements in the vector
inline double hadd_pd( __m256d vec )
{
	__m128d v2 = _mm256_extractf128_pd( vec, 1 );
	v2 = _mm_add_pd( v2, _mm256_castpd256_pd128( vec ) );
	v2 = _mm_add_sd( v2, _mm_unpackhi_pd( v2, v2 ) );
	return _mm_cvtsd_f64( v2 );
}
#endif

// Assuming class E has 2 integers, like std::pair<int,int> or similar, conditionally swap them making sure elt[0] <= elt[1]
// Only works correctly for non-negative integers
template<class E>
inline void sortInt2( E& elt )
{
	static_assert( sizeof( E ) >= 8 );
	uint64_t u = *(const uint64_t*)&elt;
	uint64_t flip = _rotr64( u, 32 );
	u = std::max( u, flip );  //< that thing compiles into a conditional move, no branches
	*(uint64_t*)&elt = u;
}

// Assuming class E has 2 integers, compare them for a < b
// Only works correctly for non-negative integers
template<class E>
inline bool compareInt2( const E& a, const E& b )
{
	static_assert( sizeof( E ) >= 8 );
	uint64_t e1 = *(const uint64_t*)&a;
	uint64_t e2 = *(const uint64_t*)&b;
	e1 = _rotr64( e1, 32 );
	e2 = _rotr64( e2, 32 );
	return e1 < e2;
}