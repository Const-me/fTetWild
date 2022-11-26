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

}

void floatTetWild::AMIPS_hessian_v2( const std::array<double, 12>& arr, Matrix3& result_0 )
{
	// Load slices of 3 elements into 4 vectors
	const __m256d v0 = _mm256_loadu_pd( &arr[ 0 ] );
	const __m256d v1 = _mm256_loadu_pd( &arr[ 3 ] );
	const __m256d v2 = _mm256_loadu_pd( &arr[ 6 ] );
	const __m256d v3 = AvxMath::loadDouble3( &arr[ 9 ] );

	const __m256d m1 = _mm256_broadcast_sd( &s_magic.m1 );
	const __m256d t00 = _mm256_mul_pd( m1, v0 );
	const __m256d t01 = _mm256_mul_pd( m1, v3 );

	const __m256d m2 = _mm256_broadcast_sd( &s_magic.m2 );
	const __m256d t02 = _mm256_mul_pd( m2, v1 );

	const __m256d m3 = _mm256_broadcast_sd( &s_magic.m3 );
	const __m256d t03 = _mm256_mul_pd( m3, v0 );
	const __m256d t04 = _mm256_mul_pd( m3, v1 );
	const __m256d t05 = _mm256_mul_pd( m3, v3 );

	const __m256d m4 = _mm256_broadcast_sd( &s_magic.m4 );
	const __m256d t06 = _mm256_mul_pd( m4, v2 );

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

	const __m256d t07 = _mm256_sub_pd( v0, v3 );
	STORE( t07 );

	const __m256d t08 = _mm256_add_pd( _mm256_sub_pd( t00, t02 ), t01 );
	STORE( t08 );

	const __m256d t09 = _mm256_add_pd( _mm256_sub_pd( _mm256_add_pd( t03, t04 ), t06 ), t05 );
	STORE( t09 );

	const __m256d cp1 = _mm256_mul_pd( permute_yxx( t08 ), permute_zzy( t09 ) );
	const __m256d cp2 = _mm256_mul_pd( permute_zzy( t08 ), permute_yxx( t09 ) );
	STORE( cp1 );
	STORE( cp2 );

	const __m256d cp = _mm256_sub_pd( cp1, cp2 );
	STORE( cp );

	const __m256d prod = _mm256_mul_pd( t07, cp );
	STORE( prod );

	const double helper_54 = prod_z + prod_x - prod_y;
	const double helper_55 = pow2( helper_54 );
	const double helper_56 = 1.0 / cubicRoot( helper_55 );
	const double helper_57 = 1.0 * v2_x - 3.0 * v0_x + 1.0 * v1_x + 1.0 * v3_x;
	const double helper_58 = magic5 * v1_y;
	const double helper_59 = magic5 * v2_y;
	const double helper_60 = helper_58 - helper_59;
	const double helper_61 = t07_z * helper_60;
	const double helper_62 = magic5 * v1_z - magic5 * v2_z;
	const double helper_63 = t07_y * helper_62;
	const double helper_64 = cp_x + helper_61 - helper_63;
	const double helper_65 = 1.33333333333333 / helper_54;
	const double helper_66 = 1.0 / helper_55;
	const double helper_67 = 0.5 * v2_x + 0.5 * v1_x;
	const double helper_68 = -1.5 * v0_x + helper_67 + 0.5 * v3_x;
	const double helper_69 = 0.5 * v0_x + helper_67 - 1.5 * v3_x;
	const double helper_70 = -1.5 * v2_x + 0.5 * v0_x + 0.5 * v1_x + 0.5 * v3_x;
	const double helper_71 = 0.5 * v2_x + 0.5 * v0_x - 1.5 * v1_x + 0.5 * v3_x;
	const double helper_72 = 0.5 * v1_y + 0.5 * v2_y;
	const double helper_73 = -1.5 * v0_y + 0.5 * v3_y + helper_72;
	const double helper_74 = 0.5 * v0_y - 1.5 * v3_y + helper_72;
	const double helper_75 = 0.5 * v0_y + 0.5 * v1_y - 1.5 * v2_y + 0.5 * v3_y;
	const double helper_76 = 0.5 * v0_y - 1.5 * v1_y + 0.5 * v2_y + 0.5 * v3_y;
	const double helper_77 = 0.5 * v1_z + 0.5 * v2_z;
	const double helper_78 = -1.5 * v0_z + 0.5 * v3_z + helper_77;
	const double helper_79 = 0.5 * v0_z - 1.5 * v3_z + helper_77;
	const double helper_80 = 0.5 * v0_z + 0.5 * v3_z + 0.5 * v1_z - 1.5 * v2_z;
	const double helper_81 = 0.5 * v0_z + 0.5 * v3_z - 1.5 * v1_z + 0.5 * v2_z;
	const double helper_82 = v0_z * helper_78 + v0_y * helper_73 + v1_y * helper_76 + v2_y * helper_75 + v3_y * helper_74 +
							 v3_z * helper_79 + v2_x * helper_70 + v1_z * helper_81 + v2_z * helper_80 + v0_x * helper_68 +
							 v1_x * helper_71 + helper_69 * v3_x;
	const double helper_83 = 0.444444444444444 * helper_66 * helper_82;
	const double helper_84 = helper_66 * helper_82;
	const double helper_85 = -prod_z - prod_x + prod_y;
	const double helper_86 = 1.0 / helper_85;
	const double helper_87 = helper_86 / cubicRoot( pow2( helper_85 ) );
	const double helper_88 = magic5 * v1_x;
	const double helper_89 = magic5 * v2_x;
	const double helper_90 = helper_88 - helper_89;
	const double helper_91 = 0.666666666666667 * t08_x * t09_z + 0.666666666666667 * t07_z * helper_90 - 0.666666666666667 * t09_x * t08_z -
							 0.666666666666667 * t07_x * helper_62;
	const double helper_92 = -3.0 * v0_y + 1.0 * v1_y + 1.0 * v2_y + 1.0 * v3_y;
	const double helper_93 = -v0_y + v3_y;
	const double helper_94 = -v0_z + v3_z;
	const double helper_95 = -t00_y + t02_y - t01_y;
	const double helper_96 = -t03_z - t04_z + t06_z - t05_z;
	const double helper_97 = -t00_z + t02_z - t01_z;
	const double helper_98 = -t03_y - t04_y + t06_y - t05_y;
	const double helper_99 = -0.666666666666667 * helper_60 * helper_94 + 0.666666666666667 * helper_62 * helper_93 +
							 0.666666666666667 * helper_95 * helper_96 - 0.666666666666667 * helper_97 * helper_98;
	const double helper_100 = t07_z * helper_90;
	const double helper_101 = t07_x * helper_62;
	const double helper_102 = helper_100 - helper_101 + cp_y;
	const double helper_103 = -helper_60 * helper_94 + helper_62 * helper_93 + helper_95 * helper_96 - helper_97 * helper_98;
	const double helper_104 = 0.444444444444444 * helper_102 * helper_103 * helper_82 * helper_86 + helper_57 * helper_91 - helper_92 * helper_99;
	const double helper_105 =
	  1.85037170770859e-17 * v0_z * helper_78 + 1.85037170770859e-17 * v0_y * helper_73 + 1.85037170770859e-17 * v1_y * helper_76 +
	  1.85037170770859e-17 * v2_y * helper_75 + 1.85037170770859e-17 * v3_y * helper_74 + 1.85037170770859e-17 * v3_z * helper_79 +
	  1.85037170770859e-17 * v2_x * helper_70 + 1.85037170770859e-17 * v1_z * helper_81 + 1.85037170770859e-17 * v2_z * helper_80 +
	  1.85037170770859e-17 * v0_x * helper_68 + 1.85037170770859e-17 * v1_x * helper_71 + 1.85037170770859e-17 * helper_69 * v3_x;
	const double helper_106 = helper_64 * helper_82 * helper_86;
	const double helper_107 = -0.666666666666667 * t08_x * t09_y + 0.666666666666667 * t08_y * t09_x +
							  0.666666666666667 * t07_x * helper_60 - 0.666666666666667 * t07_y * helper_90;
	const double helper_108 = -3.0 * v0_z + 1.0 * v3_z + 1.0 * v1_z + 1.0 * v2_z;
	const double helper_109 = -cp1_z + cp2_z + t07_x * helper_60 - t07_y * helper_90;
	const double helper_110 = 0.444444444444444 * helper_109 * helper_82 * helper_86;
	const double helper_111 = helper_103 * helper_110 + helper_107 * helper_57 - helper_108 * helper_99;
	const double helper_112 = -v0_x + v3_x;
	const double helper_113 = -helper_88 + helper_89;
	const double helper_114 = -t00_x + t02_x - t01_x;
	const double helper_115 = -t03_x - t04_x + t06_x - t05_x;
	const double helper_116 = helper_82 * helper_86 * ( helper_112 * helper_62 + helper_113 * helper_94 + helper_114 * helper_96 - helper_115 * helper_97 );
	const double helper_117 = -helper_100 + helper_101 - cp1_y + cp2_y;
	const double helper_118 = -helper_102 * helper_110 + helper_107 * helper_92 + helper_108 * helper_91;
	const double helper_119 =
	  helper_82 * helper_86 * ( helper_112 * ( -helper_58 + helper_59 ) - helper_113 * helper_93 - helper_114 * helper_98 + helper_115 * helper_95 );
	result_0( 0, 0 ) = helper_56 * ( helper_57 * helper_64 * helper_65 - pow2( helper_64 ) * helper_83 +
									 0.666666666666667 * helper_64 * helper_84 * ( -cp1_x + cp2_x - helper_61 + helper_63 ) + 3.0 );
	result_0( 0, 1 ) = helper_87 * ( helper_104 - helper_105 * v1_z + helper_106 * helper_91 );
	result_0( 0, 2 ) = helper_87 * ( helper_106 * helper_107 + helper_111 );
	result_0( 1, 0 ) = helper_87 * ( helper_104 + helper_116 * helper_99 );
	result_0( 1, 1 ) = helper_56 * ( -pow2( helper_117 ) * helper_83 + helper_117 * helper_65 * helper_92 + helper_117 * helper_84 * helper_91 + 3.0 );
	result_0( 1, 2 ) = helper_87 * ( -helper_105 * v1_x - helper_107 * helper_116 + helper_118 );
	result_0( 2, 0 ) = helper_87 * ( -helper_105 * v1_y + helper_111 + helper_119 * helper_99 );
	result_0( 2, 1 ) = helper_87 * ( helper_118 - helper_119 * helper_91 );
	result_0( 2, 2 ) = helper_56 * ( -helper_108 * helper_109 * helper_65 - 1.11111111111111 * pow2( helper_109 ) * helper_84 + 3.0 );
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