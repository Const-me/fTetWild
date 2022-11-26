#include "stdafx.h"
#include "LocalOperations.h"
#include "LocalOperations2.h"
#include <Utils/lowLevel.h>

namespace
{
	inline double pow2( double x )
	{
		return x * x;
	}

	inline double cubicRoot( double x )
	{
		return floatTetWild::cbrt( x );
	}

	constexpr double magic1 = 0.577350269189626;  // sqrt(1/3)
	constexpr double magic2 = 1.15470053837925;	  // sqrt(4/3)
	constexpr double magic3 = 0.408248290463863;  // sqrt(1/6)
	constexpr double magic4 = 1.22474487139159;	  // sqrt(3/2)
	constexpr double magic5 = 0.707106781186548;  // sqrt(0.5)

	struct sMagicNumbers
	{
		const double m1 = magic1;
		const double m2 = magic2;
		const double m3 = magic3;
		const double m4 = magic4;
		const double m5 = magic5;
		const double minus3 = -3.0;
		const double half = 0.f;
	};
	static const alignas( 64 ) sMagicNumbers s_magic;

	template<int i>
	inline double extract( __m256d v );

	template<>
	inline double extract<0>( __m256d v )
	{
		return _mm256_cvtsd_f64( v );
	}

	template<>
	inline double extract<1>( __m256d v )
	{
		__m128d low = _mm256_castpd256_pd128( v );
		low = _mm_unpackhi_pd( low, low );
		return _mm_cvtsd_f64( low );
	}

	template<>
	inline double extract<2>( __m256d v )
	{
		__m128d high = _mm256_extractf128_pd( v, 1 );
		return _mm_cvtsd_f64( high );
	}

	inline __m256d add( __m256d a, __m256d b )
	{
		return _mm256_add_pd( a, b );
	}

	inline __m256d sub( __m256d a, __m256d b )
	{
		return _mm256_sub_pd( a, b );
	}

	inline __m256d mul( __m256d a, __m256d b )
	{
		return _mm256_mul_pd( a, b );
	}

#define STORE( v )                        \
	const double v##_x = extract<0>( v ); \
	const double v##_y = extract<1>( v ); \
	const double v##_z = extract<2>( v );

#ifdef __AVX2__
	inline __m256d permute_yxx( __m256d vec )
	{
		return _mm256_permute4x64_pd( vec, _MM_SHUFFLE( 3, 0, 0, 1 ) );
	}
	inline __m256d permute_zzy( __m256d vec )
	{
		return _mm256_permute4x64_pd( vec, _MM_SHUFFLE( 3, 1, 2, 2 ) );
	}
#else
#error Not currently implemented
#endif	// __AVX2__

}  // namespace

void floatTetWild::AMIPS_hessian_v2( const std::array<double, 12>& arr, Matrix3& result_0 )
{
	// Load slices of 3 elements into 4 vectors
	const __m256d v0 = _mm256_loadu_pd( &arr[ 0 ] );
	const __m256d v1 = _mm256_loadu_pd( &arr[ 3 ] );
	const __m256d v2 = _mm256_loadu_pd( &arr[ 6 ] );
	const __m256d v3 = AvxMath::loadDouble3( &arr[ 9 ] );

	const __m256d m1 = _mm256_broadcast_sd( &s_magic.m1 );
	const __m256d t00 = mul( m1, v0 );
	const __m256d t01 = mul( m1, v3 );

	const __m256d m2 = _mm256_broadcast_sd( &s_magic.m2 );
	const __m256d t02 = mul( m2, v1 );

	const __m256d m3 = _mm256_broadcast_sd( &s_magic.m3 );
	const __m256d t03 = mul( m3, v0 );
	const __m256d t04 = mul( m3, v1 );
	const __m256d t05 = mul( m3, v3 );

	const __m256d m4 = _mm256_broadcast_sd( &s_magic.m4 );
	const __m256d t06 = mul( m4, v2 );

	STORE( v0 );
	STORE( v1 );
	STORE( v2 );
	STORE( v3 );

	STORE( t00 );
	STORE( t01 );
	STORE( t02 );
	STORE( t03 );
	STORE( t04 );
	STORE( t05 );
	STORE( t06 );

	const __m256d t07 = sub( v0, v3 );
	STORE( t07 );

	const __m256d t08 = add( sub( t00, t02 ), t01 );
	STORE( t08 );

	const __m256d t09 = add( sub( add( t03, t04 ), t06 ), t05 );
	STORE( t09 );

	const __m256d cp1 = mul( permute_yxx( t08 ), permute_zzy( t09 ) );
	const __m256d cp2 = mul( permute_zzy( t08 ), permute_yxx( t09 ) );
	STORE( cp1 );
	STORE( cp2 );

	const __m256d cp = sub( cp1, cp2 );
	STORE( cp );

	const __m256d prod = mul( t07, cp );
	STORE( prod );

	const double st0 = prod_z + prod_x - prod_y;
	const double st1 = pow2( st0 );
	const double st2 = 1.33333333333333 / st0;
	const double st3 = 1.0 / cubicRoot( st1 );
	const double st4 = 1.0 / st1;

	const double helper_57 = v2_x - 3.0 * v0_x + v1_x + v3_x;

	const __m256d m5 = _mm256_broadcast_sd( &s_magic.m5 );
	const __m256d t10 = sub( mul( m5, v1 ), mul( m5, v2 ) );
	STORE( t10 );

	const double helper_61 = t07_z * t10_y;
	const double helper_63 = t07_y * t10_z;
	const double helper_64 = cp_x + helper_61 - helper_63;

	const double helper_100 = t07_z * t10_x;
	const double helper_101 = t07_x * t10_z;
	const double helper_102 = helper_100 - helper_101 + cp_y;

	const __m256d t11 = add( v1, v2 );
	STORE( t11 );

	const double t12_x = ( -3 * v0_x + t11_x + v3_x ) * 0.5;
	const double t12_y = ( -3 * v0_y + v3_y + t11_y ) * 0.5;
	const double t12_z = ( -3 * v0_z + v3_z + t11_z ) * 0.5;

	const double t13_x = ( v0_x + t11_x - 3 * v3_x ) * 0.5;
	const double t13_y = ( v0_y - 3 * v3_y + t11_y ) * 0.5;
	const double t13_z = ( v0_z - 3 * v3_z + t11_z ) * 0.5;

	const double t14_x = ( -3 * v2_x + v0_x + v1_x + v3_x ) * 0.5;
	const double t14_y = ( v0_y + v1_y - 3 * v2_y + v3_y ) * 0.5;
	const double t14_z = ( v0_z + v3_z + v1_z - 3 * v2_z ) * 0.5;

	const double t15_x = ( v2_x + v0_x - 3 * v1_x + v3_x ) * 0.5;
	const double t15_y = ( v0_y - 3 * v1_y + v2_y + v3_y ) * 0.5;
	const double t15_z = ( v0_z + v3_z - 3 * v1_z + v2_z ) * 0.5;

	const double product1 = v0_z * t12_z + v0_y * t12_y + v1_y * t15_y + v2_y * t14_y + v3_y * t13_y + v3_z * t13_z + v2_x * t14_x + v1_z * t15_z +
							v2_z * t14_z + v0_x * t12_x + v1_x * t15_x + t13_x * v3_x;

	const double helper_83 = 0.444444444444444 * st4 * product1;
	const double helper_84 = st4 * product1;
	const double helper_85 = -st0;
	const double helper_86 = 1.0 / helper_85;
	const double helper_87 = helper_86 / cubicRoot( pow2( helper_85 ) );
	const double helper_91 =
	  0.666666666666667 * t08_x * t09_z + 0.666666666666667 * t07_z * t10_x - 0.666666666666667 * t09_x * t08_z - 0.666666666666667 * t07_x * t10_z;
	const double helper_92 = -3.0 * v0_y + 1.0 * v1_y + 1.0 * v2_y + 1.0 * v3_y;

	const double t16_x = v3_x - v0_x;
	const double t16_y = v3_y - v0_y;
	const double t16_z = v3_z - v0_z;

	const double t17_x = -t00_x + t02_x - t01_x;
	const double t17_y = -t00_y + t02_y - t01_y;
	const double t17_z = -t00_z + t02_z - t01_z;

	const double t18_x = -t03_x - t04_x + t06_x - t05_x;
	const double t18_y = -t03_y - t04_y + t06_y - t05_y;
	const double t18_z = -t03_z - t04_z + t06_z - t05_z;

	const double helper_99 =
	  -0.666666666666667 * t10_y * t16_z + 0.666666666666667 * t10_z * t16_y + 0.666666666666667 * t17_y * t18_z - 0.666666666666667 * t17_z * t18_y;
	const double helper_103 = -t10_y * t16_z + t10_z * t16_y + t17_y * t18_z - t17_z * t18_y;
	const double helper_104 = 0.444444444444444 * helper_102 * helper_103 * product1 * helper_86 + helper_57 * helper_91 - helper_92 * helper_99;
	const double helper_105 = 1.85037170770859e-17 * v0_z * t12_z + 1.85037170770859e-17 * v0_y * t12_y + 1.85037170770859e-17 * v1_y * t15_y +
							  1.85037170770859e-17 * v2_y * t14_y + 1.85037170770859e-17 * v3_y * t13_y + 1.85037170770859e-17 * v3_z * t13_z +
							  1.85037170770859e-17 * v2_x * t14_x + 1.85037170770859e-17 * v1_z * t15_z + 1.85037170770859e-17 * v2_z * t14_z +
							  1.85037170770859e-17 * v0_x * t12_x + 1.85037170770859e-17 * v1_x * t15_x + 1.85037170770859e-17 * t13_x * v3_x;

	const double helper_106 = helper_64 * product1 * helper_86;
	const double helper_107 =
	  -0.666666666666667 * t08_x * t09_y + 0.666666666666667 * t08_y * t09_x + 0.666666666666667 * t07_x * t10_y - 0.666666666666667 * t07_y * t10_x;
	const double helper_108 = -3.0 * v0_z + 1.0 * v3_z + 1.0 * v1_z + 1.0 * v2_z;
	const double helper_109 = -cp_z + t07_x * t10_y - t07_y * t10_x;
	const double helper_110 = 0.444444444444444 * helper_109 * product1 * helper_86;
	const double helper_111 = helper_103 * helper_110 + helper_107 * helper_57 - helper_108 * helper_99;
	const double helper_116 = product1 * helper_86 * ( t16_x * t10_z - t10_x * t16_z + t17_x * t18_z - t18_x * t17_z );
	const double helper_117 = -helper_100 + helper_101 - cp1_y + cp2_y;
	const double helper_118 = -helper_102 * helper_110 + helper_107 * helper_92 + helper_108 * helper_91;
	const double helper_119 = product1 * helper_86 * ( -t16_x * t10_y + t10_x * t16_y - t17_x * t18_y + t18_x * t17_y );
	result_0( 0, 0 ) = st3 * ( helper_57 * helper_64 * st2 - pow2( helper_64 ) * helper_83 +
							   0.666666666666667 * helper_64 * helper_84 * ( -cp_x - helper_61 + helper_63 ) + 3.0 );
	result_0( 0, 1 ) = helper_87 * ( helper_104 - helper_105 * v1_z + helper_106 * helper_91 );
	result_0( 0, 2 ) = helper_87 * ( helper_106 * helper_107 + helper_111 );
	result_0( 1, 0 ) = helper_87 * ( helper_104 + helper_116 * helper_99 );
	result_0( 1, 1 ) = st3 * ( -pow2( helper_117 ) * helper_83 + helper_117 * st2 * helper_92 + helper_117 * helper_84 * helper_91 + 3.0 );
	result_0( 1, 2 ) = helper_87 * ( -helper_105 * v1_x - helper_107 * helper_116 + helper_118 );
	result_0( 2, 0 ) = helper_87 * ( -helper_105 * v1_y + helper_111 + helper_119 * helper_99 );
	result_0( 2, 1 ) = helper_87 * ( helper_118 - helper_119 * helper_91 );
	result_0( 2, 2 ) = st3 * ( -helper_108 * helper_109 * st2 - 1.11111111111111 * pow2( helper_109 ) * helper_84 + 3.0 );
}

void floatTetWild::AMIPS_hessian( const std::array<Scalar, 12>& T, Matrix3& result_0 )
{
#if 0
	AMIPS_hessian_v2( T, result_0 );
#else
	Matrix3 matOld, matNew;
	AMIPS_hessian_v1( T, matOld );
	AMIPS_hessian_v2( T, matNew );
	Matrix3 diff = matNew - matOld;
	__debugbreak();
	result_0 = matOld;
#endif
}