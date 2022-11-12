#pragma once
#ifdef __AVX__
#include "AvxMath.h"
#else
#include <emmintrin.h>
#include <smmintrin.h>
#endif

namespace Simd
{
	struct alignas( 32 ) Vec3
	{
#ifdef __AVX__
		__m256d v;

		Vec3 operator-( const Vec3& b ) const
		{
			return Vec3 { _mm256_sub_pd( v, b ) };
		}
		Vec3 operator*( double s ) const
		{
			__m256d sv = _mm256_set1_pd( s );
			return Vec3 { _mm256_mul_pd( v, sv ) };
		}
		operator __m256d() const
		{
			return v;
		}
#else
		__m128d xy, z;

		Vec3 operator-( const Vec3& b ) const
		{
			return Vec3 { _mm_sub_pd( xy, b.xy ), _mm_sub_pd( z, b.z ) };
		}

		Vec3 operator*( double s ) const
		{
			__m128d sv = _mm_set1_pd( s );
			return Vec3 { _mm_mul_pd( xy, sv ), _mm_mul_pd( z, sv ) };
		}

		inline __m128d yz() const
		{
			return _mm_shuffle_pd( xy, z, _MM_SHUFFLE2( 0, 1 ) );
		}
		inline __m128d zx() const
		{
			return _mm_unpacklo_pd( z, xy );
		}
		inline __m128d yy() const
		{
			return _mm_unpackhi_pd( xy, xy );
		}
#endif
	};

#ifdef __AVX__
	inline Vec3 load3( const double* rsi )
	{
		return Vec3 { AvxMath::loadDouble3( rsi ) };
	}
	inline void store3( double* rdi, Vec3 v )
	{
		AvxMath::storeDouble3( rdi, v );
	}

	inline double dot( const Vec3& a, const Vec3& b )
	{
		return AvxMath::vector3DotScalar( a, b );
	}

	inline Vec3 cross( const Vec3& a, const Vec3& b )
	{
		return Vec3 { AvxMath::vector3Cross( a, b ) };
	}
#else
	inline Vec3 load3( const double* rsi )
	{
		return Vec3 { _mm_loadu_pd( rsi ), _mm_load_sd( rsi + 2 ) };
	}

	inline void store3( double* rdi, Vec3 v )
	{
		_mm_storeu_pd( rdi, v.xy );
		_mm_store_sd( rdi + 2, v.z );
	}

	inline double dot( const Vec3& a, const Vec3& b )
	{
		__m128d r = _mm_dp_pd( a.xy, b.xy, 0b00110001 );
		__m128d zz = _mm_mul_sd( a.z, b.z );
		r = _mm_add_sd( r, zz );
		return _mm_cvtsd_f64( r );
	}

	inline Vec3 cross( const Vec3& a, const Vec3& b )
	{
		Vec3 res;
		// a.yzx * b.zxy - a.zxy - b*yzx
		res.xy = _mm_sub_pd( _mm_mul_pd( a.yz(), b.zx() ), _mm_mul_pd( a.zx(), b.yz() ) );
		res.z = _mm_sub_pd( _mm_mul_pd( a.xy, b.yy() ), _mm_mul_pd( a.yy(), b.xy ) );
		return res;
	}
#endif
}  // namespace Simd