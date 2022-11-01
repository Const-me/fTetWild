#pragma once
#include <immintrin.h>

// Load 12 bytes from memory
inline __m128i loadInt3( const void* rsi )
{
	__m128i v = _mm_loadu_si64( rsi );
	v = _mm_insert_epi32( v, ( (const int*)( rsi ) )[ 2 ], 2 );
	return v;
}

// Make a mask for _mm_insert_ps instruction
constexpr int insertMask( int source, int dest, int zero = 0 )
{
	assert( source >= 0 && source < 4 );
	assert( dest >= 0 && dest < 4 );
	assert( zero >= 0 && zero < 16 );
	return ( source << 6 ) | ( dest << 4 ) | zero;
}

// Load 3 floats from memory
inline __m128 loadFloat3( const void* rsi )
{
	__m128d d = _mm_load_sd( (const double*)rsi );
	__m128 f = _mm_castpd_ps( d );
	__m128 z = _mm_load_ss( ( (const float*)rsi ) + 2 );
	return _mm_insert_ps( f, z, insertMask( 0, 2, 0 ) );
}

inline bool vectorEqual( __m128 a, __m128 b )
{
	// compare for a != p but  we want TRUE for NAN, that's why the weird predicate
	__m128 cmp = _mm_cmp_ps( a, b, _CMP_NEQ_UQ );
	return (bool)_mm_testz_ps( cmp, cmp );
}

// Store 12 floats from XYZ lanes of 4 vectors, with 3 store instructions.
// The function writes 12 floats = 48 bytes to the pointer.
inline void storeFloat3_x4( __m128 a, __m128 b, __m128 c, __m128 d, void* dest )
{
	float* const rdi = (float*)dest;

	// a.xyz, b.x
	__m128 tmp = _mm_insert_ps( a, b, insertMask( 0, 3, 0 ) );
	_mm_storeu_ps( rdi, tmp );

	// b.yz, c.xy
	tmp = _mm_shuffle_ps( b, c, _MM_SHUFFLE( 1, 0, 2, 1 ) );
	_mm_storeu_ps( rdi + 4, tmp );

	// c.z, d.xyz
	tmp = _mm_permute_ps( d, _MM_SHUFFLE( 2, 1, 0, 0 ) );  // Permute into d.xxyz
	tmp = _mm_insert_ps( tmp, c, insertMask( 2, 0, 0 ) );
	_mm_storeu_ps( rdi + 8, tmp );
}