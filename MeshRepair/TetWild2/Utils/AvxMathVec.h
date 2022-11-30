#pragma once
#include "AvxMath.h"

namespace AvxMath
{
	// A wrapper around __m256d implementing arithmetic operators
	// The only reason for a separate structure - neither gcc nor clang support operators outside of classes.
	// For VC++ it's possible to define these operators directly on the __m256d type
	struct Vec
	{
		__m256d vec;

		Vec( __m256d v )
			: vec( v )
		{
		}

		operator __m256d() const
		{
			return vec;
		}

		Vec operator+( __m256d v ) const
		{
			return _mm256_add_pd( vec, v );
		}
		Vec operator-( __m256d v ) const
		{
			return _mm256_sub_pd( vec, v );
		}
		Vec operator*( __m256d v ) const
		{
			return _mm256_mul_pd( vec, v );
		}
		Vec operator/( __m256d v ) const
		{
			return _mm256_div_pd( vec, v );
		}
		void operator+=( __m256d v )
		{
			vec = _mm256_add_pd( vec, v );
		}
		void operator-=( __m256d v )
		{
			vec = _mm256_sub_pd( vec, v );
		}
		void operator*=( __m256d v )
		{
			vec = _mm256_mul_pd( vec, v );
		}
		void operator/=( __m256d v )
		{
			vec = _mm256_div_pd( vec, v );
		}
		Vec operator*( double s ) const
		{
			return _mm256_mul_pd( vec, _mm256_set1_pd( s ) );
		}
		void operator*=( double s )
		{
			vec = _mm256_mul_pd( vec, _mm256_set1_pd( s ) );
		}
		Vec operator/( double s ) const
		{
			return _mm256_div_pd( vec, _mm256_set1_pd( s ) );
		}
		void operator/=( double s )
		{
			vec = _mm256_div_pd( vec, _mm256_set1_pd( s ) );
		}
	};

	inline Vec normalize( __m256d vec )
	{
		return vector3Normalize( vec );
	}

	inline double distance( __m256d a, __m256d b )
	{
		return vector3Length( _mm256_sub_pd( a, b ) );
	}

	inline double dot( __m256d a, __m256d b )
	{
		return vector3DotScalar( a, b );
	}

	inline Vec cross( __m256d a, __m256d b )
	{
		return vector3Cross( a, b );
	}
}  // namespace AvxMath