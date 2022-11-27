#include "stdafx.h"
#include "LocalOperations.h"
#include "LocalOperations2.h"
#include <Utils/lowLevel.h>
#include "LocalOperationsUtils.h"

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
		const double half = 0.5;
	};
	static const alignas( 64 ) sMagicNumbers s_magic;

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
	const double root = cubicRoot( st1 );
	const double st3 = 1.0 / root;
	const double st4 = 1.0 / st1;
	const double st5 = -1.0 / st0;
	const double st6 = st5 / root;

	const __m256d m5 = _mm256_broadcast_sd( &s_magic.m5 );
	const __m256d t10 = sub( mul( m5, v1 ), mul( m5, v2 ) );
	STORE( t10 );

	const double helper_100 = t07_z * t10_x;
	const double helper_101 = t07_x * t10_z;
	const double helper_102 = helper_100 - helper_101 + cp_y;

	const double t19_x = -cp_x - t07_z * t10_y + t07_y * t10_z;
	// cp = cp1 - cp2, but simplifying would break binary equality due to different summation order
	const double t19_y = -helper_100 + helper_101 - cp1_y + cp2_y;
	const double t19_z = -cp_z + t07_x * t10_y - t07_y * t10_x;

	const __m256d t11 = add( v1, v2 );
	STORE( t11 );

	const double t12_x = -3 * v0_x + t11_x + v3_x;
	const double t12_y = -3 * v0_y + v3_y + t11_y;
	const double t12_z = -3 * v0_z + v3_z + t11_z;

	const double t13_x = v0_x + t11_x - 3 * v3_x;
	const double t13_y = v0_y - 3 * v3_y + t11_y;
	const double t13_z = v0_z - 3 * v3_z + t11_z;

	const double t14_x = -3 * v2_x + v0_x + v1_x + v3_x;
	const double t14_y = v0_y + v1_y - 3 * v2_y + v3_y;
	const double t14_z = v0_z + v3_z + v1_z - 3 * v2_z;

	const double t15_x = v2_x + v0_x - 3 * v1_x + v3_x;
	const double t15_y = v0_y - 3 * v1_y + v2_y + v3_y;
	const double t15_z = v0_z + v3_z - 3 * v1_z + v2_z;

	const double product1 = ( v0_z * t12_z + v0_y * t12_y + v1_y * t15_y + v2_y * t14_y + v3_y * t13_y + v3_z * t13_z + v2_x * t14_x + v1_z * t15_z +
							  v2_z * t14_z + v0_x * t12_x + v1_x * t15_x + t13_x * v3_x ) *
							0.5;

	const double st7 = st4 * product1;

	const double t22_x = -3.0 * v0_x + v2_x + v1_x + v3_x;
	const double t22_y = -3.0 * v0_y + v1_y + v2_y + v3_y;
	const double t23_y = -3.0 * v0_z + v3_z + v1_z + v2_z;

	const __m256d t16 = sub( v3, v0 );
	STORE( t16 );

	const __m256d t17 = sub( sub( t02, t00 ), t01 );
	STORE( t17 );

	const __m256d t18 = AvxMath::vectorNegate( add( sub( add( t03, t04 ), t06 ), t05 ) );
	STORE( t18 );

	const double t20_x =
	  -0.666666666666667 * t10_y * t16_z + 0.666666666666667 * t10_z * t16_y + 0.666666666666667 * t17_y * t18_z - 0.666666666666667 * t17_z * t18_y;
	const double t20_y =
	  0.666666666666667 * t08_x * t09_z + 0.666666666666667 * t07_z * t10_x - 0.666666666666667 * t09_x * t08_z - 0.666666666666667 * t07_x * t10_z;
	const double t20_z =
	  -0.666666666666667 * t08_x * t09_y + 0.666666666666667 * t08_y * t09_x + 0.666666666666667 * t07_x * t10_y - 0.666666666666667 * t07_y * t10_x;

	const __m256d t21_0 = mul( permute_yxx( t10 ), permute_zzy( t16 ) );
	const __m256d t21_1 = mul( permute_yxx( t16 ), permute_zzy( t10 ) );
	const __m256d t21_2 = mul( permute_yxx( t17 ), permute_zzy( t18 ) );
	const __m256d t21_3 = mul( permute_yxx( t18 ), permute_zzy( t17 ) );
	const __m256d t21 = sub( add( sub( t21_1, t21_0 ), t21_2 ), t21_3 );
	STORE( t21 );

	const double helper_104 = 0.444444444444444 * helper_102 * t21_x * product1 * st5 + t22_x * t20_y - t22_y * t20_x;

	const double helper_105 = ( 1.85037170770859e-17 * v0_z * t12_z + 1.85037170770859e-17 * v0_y * t12_y + 1.85037170770859e-17 * v1_y * t15_y +
								1.85037170770859e-17 * v2_y * t14_y + 1.85037170770859e-17 * v3_y * t13_y + 1.85037170770859e-17 * v3_z * t13_z +
								1.85037170770859e-17 * v2_x * t14_x + 1.85037170770859e-17 * v1_z * t15_z + 1.85037170770859e-17 * v2_z * t14_z +
								1.85037170770859e-17 * v0_x * t12_x + 1.85037170770859e-17 * v1_x * t15_x + 1.85037170770859e-17 * t13_x * v3_x ) *
							  0.5;

	const double helper_106 = -t19_x * product1 * st5;
	const double helper_83 = 0.444444444444444 * st4 * product1;
	const double helper_110 = 0.444444444444444 * t19_z * product1 * st5;
	const double helper_111 = t21_x * helper_110 + t20_z * t22_x - t23_y * t20_x;
	const double helper_116 = product1 * st5 * t21_y;
	const double helper_118 = -helper_102 * helper_110 + t20_z * t22_y + t23_y * t20_y;
	const double helper_119 = product1 * st5 * t21_z;

	const double t24_x = -pow2( t19_x );
	const double t24_y = -pow2( t19_y );
	const double t24_z = -pow2( t19_z );

	const double diag_x = t22_x * -t19_x * st2 + t24_x * helper_83 - 0.666666666666667 * t19_x * st7 * t19_x + 3.0;
	const double diag_y = t24_y * helper_83 + t19_y * st2 * t22_y + t19_y * st7 * t20_y + 3.0;
	const double diag_z = -t23_y * t19_z * st2 + 1.11111111111111 * t24_z * st7 + 3.0;

	result_0( 0, 0 ) = st3 * diag_x;
	result_0( 0, 1 ) = st6 * ( helper_104 - helper_105 * v1_z + helper_106 * t20_y );
	result_0( 0, 2 ) = st6 * ( helper_106 * t20_z + helper_111 );
	result_0( 1, 0 ) = st6 * ( helper_104 + helper_116 * t20_x );
	result_0( 1, 1 ) = st3 * diag_y;
	result_0( 1, 2 ) = st6 * ( -helper_105 * v1_x - t20_z * helper_116 + helper_118 );
	result_0( 2, 0 ) = st6 * ( -helper_105 * v1_y + helper_111 - helper_119 * t20_x );
	result_0( 2, 1 ) = st6 * ( helper_118 + helper_119 * t20_y );
	result_0( 2, 2 ) = st3 * diag_z;
}

void floatTetWild::AMIPS_hessian_v3( const std::array<double, 12>& arr, Matrix3& result_0 )
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

	STORE( v1 );

	const __m256d t07 = sub( v0, v3 );
	STORE( t07 );

	const __m256d t08 = add( sub( t00, t02 ), t01 );
	STORE( t08 );

	const __m256d t09 = add( sub( add( t03, t04 ), t06 ), t05 );
	STORE( t09 );

	const __m256d cp = customProduct( t08, t09 );

	const __m256d prod = mul( t07, cp );

	const double st0 = horizontalAddSub( prod );
	const double st1 = pow2( st0 );
	const double st2 = 1.33333333333333 / st0;
	const double root = cubicRoot( st1 );
	const double st3 = 1.0 / root;
	const double st4 = 1.0 / st1;
	const double st5 = -1.0 / st0;
	const double st6 = st5 / root;

	const __m256d m5 = _mm256_broadcast_sd( &s_magic.m5 );
	const __m256d t10 = mul( m5, sub( v1, v2 ) );
	STORE( t10 );

	const __m256d t11 = add( v1, v2 );

	const __m256d three = _mm256_set1_pd( 3 );
	// v3 + t11 - 3.0 * v0
	const __m256d t12 = sub( add( t11, v3 ), mul( three, v0 ) );

	// v0 + t11 - 3.0 * v3
	const __m256d t13 = sub( add( t11, v0 ), mul( three, v3 ) );

	// v0 + v1 + v3 - 3.0 * v2
	const __m256d t14 = sub( add( add( v0, v1 ), v3 ), mul( three, v2 ) );

	// v0 + v2 + v3 - 3.0 * v1
	const __m256d t15 = sub( add( add( v0, v2 ), v3 ), mul( three, v1 ) );

	const double product1 = hadd12( mul( v0, t12 ), mul( v1, t15 ), mul( v2, t14 ), mul( v3, t13 ) ) * 0.5;

	const double st7 = st4 * product1;

	STORE( t12 );

	const __m256d t16 = sub( v3, v0 );
	STORE( t16 );

	const __m256d t17 = sub( sub( t02, t00 ), t01 );
	STORE( t17 );

	const __m256d t18 = AvxMath::vectorNegate( add( sub( add( t03, t04 ), t06 ), t05 ) );
	STORE( t18 );

	const double t20_x = 0.666666666666667 * ( -t10_y * t16_z + t10_z * t16_y + t17_y * t18_z - t17_z * t18_y );
	const double t20_y = 0.666666666666667 * ( t08_x * t09_z - t08_z * t09_x + t07_z * t10_x - t07_x * t10_z );
	const double t20_z = 0.666666666666667 * ( -t08_x * t09_y + t08_y * t09_x + t07_x * t10_y - t07_y * t10_x );

	const __m256d t21 = add( customProduct( t16, t10 ), customProduct( t17, t18 ) );
	STORE( t21 );

	const double helper_104 = 0.444444444444444 * t21_y * t21_x * product1 * st5 + t12_x * t20_y - t12_y * t20_x;

	const double helper_105 = 1.85037170770859e-17 * product1;

	const double st8 = product1 * st5;
	const double t25_x = t21_x * st8;
	const double t25_y = t21_y * st8;
	const double t25_z = t21_z * st8;

	const double helper_83 = 0.444444444444444 * st4 * product1;
	const double helper_110 = -0.444444444444444 * t21_z * product1 * st5;
	const double helper_111 = t21_x * helper_110 + t20_z * t12_x - t12_z * t20_x;
	const double helper_118 = -t21_y * helper_110 + t20_z * t12_y + t12_z * t20_y;

	const __m256d t24 = mul( t21, t21 );
	STORE( t24 );

	const double diag_x = -t24_x * ( helper_83 + 0.666666666666667 * st7 ) + t12_x * t21_x * st2 + 3.0;
	const double diag_y = -t24_y * helper_83 - t12_y * t21_y * st2 - t21_y * st7 * t20_y + 3.0;
	const double diag_z = -t24_z * ( 1.11111111111111 * st7 ) + t12_z * t21_z * st2 + 3.0;

	result_0( 0, 0 ) = st3 * diag_x;
	result_0( 0, 1 ) = st6 * ( helper_104 - helper_105 * v1_z + t25_x * t20_y );
	result_0( 0, 2 ) = st6 * ( t25_x * t20_z + helper_111 );
	result_0( 1, 0 ) = st6 * ( helper_104 + t25_y * t20_x );
	result_0( 1, 1 ) = st3 * diag_y;
	result_0( 1, 2 ) = st6 * ( -helper_105 * v1_x - t20_z * t25_y + helper_118 );
	result_0( 2, 0 ) = st6 * ( -helper_105 * v1_y + helper_111 - t25_z * t20_x );
	result_0( 2, 1 ) = st6 * ( helper_118 + t25_z * t20_y );
	result_0( 2, 2 ) = st3 * diag_z;
}

void floatTetWild::AMIPS_hessian_v4( const std::array<double, 12>& arr, Matrix3& result_0 )
{
	// Load slices of 3 elements into 4 vectors
	Vec v0 = _mm256_loadu_pd( &arr[ 0 ] );
	Vec v1 = _mm256_loadu_pd( &arr[ 3 ] );
	Vec v2 = _mm256_loadu_pd( &arr[ 6 ] );
	Vec v3 = AvxMath::loadDouble3( &arr[ 9 ] );
	translateToBarycenter( v0.vec, v1.vec, v2.vec, v3.vec, _mm256_set1_pd( 0.25 ) );

	const Vec m1 = broadcast( s_magic.m1 );
	const Vec t00 = m1 * v0;
	const Vec t01 = m1 * v3;

	const Vec t02 = v1 * broadcast( s_magic.m2 );

	const Vec m3 = broadcast( s_magic.m3 );
	const Vec t03 = m3 * v0;
	const Vec t04 = m3 * v1;
	const Vec t05 = m3 * v3;

	const Vec t06 = v2 * broadcast( s_magic.m4 );
	STORE( v1 );

	const Vec t07 = v0 - v3;
	STORE( t07 );

	const Vec t08 = t00 - t02 + t01;
	STORE( t08 );

	const Vec t09 = t03 + t04 - t06 + t05;
	STORE( t09 );

	using namespace AvxMath;
	const Vec cp = vector3Cross( t08, t09 );
	const Vec prod = t07 * cp;

	const double st0 = vector3HorizontalSum( prod );
	const double st1 = pow2( st0 );
	const double st2 = 1.33333333333333 / st0;
	const double root = cubicRoot( st1 );
	const double st3 = 1.0 / root;
	const double st4 = 1.0 / st1;
	const double st5 = -1.0 / st0;
	const double st6 = st5 / root;

	const Vec t10 = ( v1 - v2 ) * broadcast( s_magic.m5 );
	STORE( t10 );

	const double product1 = hadd12( v0 * v0, v1 * v1, v2 * v2, v3 * v3 ) * -2.0;
	const double st7 = st4 * product1;

	const Vec t12 = v0 * _mm256_set1_pd( -4.0 );
	STORE( t12 );

	const Vec t21 = add( vector3Cross( t08, t09 ), vector3Cross( t10, t07 ) );
	STORE( t21 );

	const double scaledProduct = 1.85037170770859e-17 * product1;

	const double st8 = product1 * st5;
	const Vec t25 = t21 * _mm256_set1_pd( st8 );
	STORE( t25 );

	const double t26_x = t21_y * t21_z * ( -4.0 / 9 ) * product1 * st5 - t21_z * t12_y * ( 6.0 / 9 ) - t12_z * t21_y * ( 6.0 / 9 );
	const double t26_y = t21_x * t21_z * ( -4.0 / 9 ) * product1 * st5 - t21_z * t12_x * ( 6.0 / 9 ) - t12_z * t21_x * ( 6.0 / 9 );
	const double t26_z = t21_y * t21_x * ( -4.0 / 9 ) * product1 * st5 - t12_x * t21_y * ( 6.0 / 9 ) - t12_y * t21_x * ( 6.0 / 9 );

	const Vec t24 = t21 * t21;
	STORE( t24 );

	const double diag_x = -t24_x * ( ( 10.0 / 9 ) * st7 ) + t12_x * t21_x * st2 + 3.0;
	const double diag_y = -t24_y * ( ( 10.0 / 9 ) * st7 ) + t12_y * t21_y * st2 + 3.0;
	const double diag_z = -t24_z * ( ( 10.0 / 9 ) * st7 ) + t12_z * t21_z * st2 + 3.0;

	result_0( 0, 0 ) = st3 * diag_x;
	result_0( 0, 1 ) = st6 * ( t26_z - scaledProduct * v1_z - t25_x * t21_y * ( 6.0 / 9 ) );
	result_0( 0, 2 ) = st6 * ( t26_y - t25_x * t21_z * ( 6.0 / 9 ) );
	result_0( 1, 0 ) = st6 * ( t26_z - t25_y * t21_x * ( 6.0 / 9 ) );
	result_0( 1, 1 ) = st3 * diag_y;
	result_0( 1, 2 ) = st6 * ( -scaledProduct * v1_x - t21_z * t25_y * ( 6.0 / 9 ) + t26_x );
	result_0( 2, 0 ) = st6 * ( -scaledProduct * v1_y + t26_y - t25_z * t21_x * ( 6.0 / 9 ) );
	result_0( 2, 1 ) = st6 * ( t26_x - t25_z * t21_y * ( 6.0 / 9 ) );
	result_0( 2, 2 ) = st3 * diag_z;
}

void floatTetWild::AMIPS_hessian( const std::array<Scalar, 12>& T, Matrix3& result_0 )
{
#if 0
	AMIPS_hessian_v4( T, result_0 );
#else
	Matrix3 matOld, matNew;
	AMIPS_hessian_v1( T, matOld );
	AMIPS_hessian_v4( T, matNew );
	Matrix3 diff = matNew - matOld;
	result_0 = matOld;
#endif
}