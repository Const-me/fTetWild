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
}  // namespace

void floatTetWild::AMIPS_jacobian_v2( const std::array<Scalar, 12>& arr, Vector3& result_0 )
{
	const double v0_x = arr[ 0 ];
	const double v0_y = arr[ 1 ];
	const double v0_z = arr[ 2 ];
	const double v1_x = arr[ 3 ];
	const double v1_y = arr[ 4 ];
	const double v1_z = arr[ 5 ];
	const double v2_x = arr[ 6 ];
	const double v2_y = arr[ 7 ];
	const double v2_z = arr[ 8 ];
	const double v3_x = arr[ 9 ];
	const double v3_y = arr[ 10 ];
	const double v3_z = arr[ 11 ];

	const double t0_x = v0_x - v3_x;
	const double t0_y = v0_y - v3_y;
	const double t0_z = v0_z - v3_z;

	const double t1_x = magic1 * v0_x - magic2 * v1_x + magic1 * v3_x;
	const double t1_y = magic1 * v0_y - magic2 * v1_y + magic1 * v3_y;
	const double t1_z = magic1 * v0_z - magic2 * v1_z + magic1 * v3_z;

	// Note different summation order of lanes, impossible to vectorize efficiently
	const double t2_x = magic3 * v0_x - magic4 * v2_x + magic3 * v1_x + magic3 * v3_x;
	const double t2_y = magic3 * v0_y + magic3 * v1_y - magic4 * v2_y + magic3 * v3_y;
	const double t2_z = magic3 * v1_z - magic4 * v2_z + magic3 * v3_z + magic3 * v0_z;

	const double prod_x = t1_y * t2_z - t1_z * t2_y;
	const double prod_y = t1_z * t2_x - t1_x * t2_z;
	const double helper_32 = t1_x * t2_y;
	const double helper_37 = t1_y * t2_x;
	const double prod_z = helper_32 - helper_37;
	const double s0 = -( prod_y * t0_y + t0_z * prod_z + t0_x * prod_x );
	const double s1 = 1.0 / cubicRoot( pow2( s0 ) );
	const double t3_x = magic5 * v1_x - magic5 * v2_x;
	const double t3_y = magic5 * v1_y - magic5 * v2_y;
	const double t3_z = magic5 * v1_z - magic5 * v2_z;
	const double t4_x = v1_x + v2_x;
	const double t4_y = v1_y + v2_y;
	const double t4_z = v1_z + v2_z;
	const double s2 = ( 0.666666666666667 * 0.5 ) *
					  ( v0_y * ( -3 * v0_y + v3_y + t4_y ) + v1_z * ( -3 * v1_z + v2_z + v3_z + v0_z ) + v2_z * ( v1_z - 3 * v2_z + v3_z + v0_z ) +
						v3_z * ( -3 * v3_z + t4_z + v0_z ) + v3_y * ( v0_y - 3 * v3_y + t4_y ) + v2_x * ( -3 * v2_x + v0_x + v1_x + v3_x ) +
						v1_y * ( v0_y + v3_y - 3 * v1_y + v2_y ) + v2_y * ( v0_y + v3_y + v1_y - 3 * v2_y ) + v0_x * ( -3 * v0_x + t4_x + v3_x ) +
						v1_x * ( v2_x + v0_x - 3 * v1_x + v3_x ) + v3_x * ( v0_x + t4_x - 3 * v3_x ) + v0_z * ( v3_z + t4_z - 3 * v0_z ) ) /
					  s0;
	result_0[ 0 ] = s1 * ( -v2_x + 3.0 * v0_x + s2 * ( t3_z * t0_y - t3_y * t0_z + t1_z * t2_y - t1_y * t2_z ) - v1_x - v3_x );
	result_0[ 1 ] = s1 * ( 3.0 * v0_y - v3_y - v1_y - v2_y + s2 * ( -prod_y + t0_z * t3_x - t0_x * t3_z ) );
	result_0[ 2 ] = s1 * ( -v1_z - v2_z - v3_z + s2 * ( -t0_y * t3_x - helper_32 + helper_37 + t0_x * t3_y ) + 3.0 * v0_z );
}

namespace
{
	struct sMagicNumbersV4
	{
		const double m1 = magic1;
		const double m2 = magic2;
		const double m3 = magic3;
		const double m4 = magic4;
		const double m5 = magic5;
		const double minus3 = -3.0;
	};
	static const alignas( 64 ) sMagicNumbersV4 s_magicv4;

	class Vec
	{
		__m256d vec;

	  public:
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
			return Vec { _mm256_add_pd( vec, v ) };
		}
		Vec operator-( __m256d v ) const
		{
			return Vec { _mm256_sub_pd( vec, v ) };
		}
		Vec operator*( __m256d v ) const
		{
			return Vec { _mm256_mul_pd( vec, v ) };
		}
	};
	inline __m256d broadcast( const double& r )
	{
		return _mm256_broadcast_sd( &r );
	}
}  // namespace

void floatTetWild::AMIPS_jacobian_v3( const std::array<Scalar, 12>& arr, Vector3& result_0 )
{
	// Load slices of 3 elements into 4 vectors
	const Vec v0 = _mm256_loadu_pd( &arr[ 0 ] );
	const Vec v1 = _mm256_loadu_pd( &arr[ 3 ] );
	const Vec v2 = _mm256_loadu_pd( &arr[ 6 ] );
	const Vec v3 = AvxMath::loadDouble3( &arr[ 9 ] );

	const Vec t0 = v0 - v3;
	const Vec magic1 = broadcast( s_magicv4.m1 );
	const Vec magic2 = broadcast( s_magicv4.m2 );
	const Vec t1 = magic1 * v0 - magic2 * v1 + magic1 * v3;

	const Vec magic3 = broadcast( s_magicv4.m3 );
	const Vec magic4 = broadcast( s_magicv4.m4 );
	const Vec t2 = magic3 * v0 + magic3 * v1 - magic4 * v2 + magic3 * v3;

	using namespace AvxMath;
	const Vec prod = vector3Cross( t1, t2 );
	const double s0 = -vector3DotScalar( prod, t0 );
	const double s1 = 1.0 / cubicRoot( pow2( s0 ) );

	const Vec t3 = ( v1 - v2 ) * broadcast( s_magicv4.m5 );
	const Vec t4 = v1 + v2;

	const Vec nm3 = broadcast( s_magicv4.minus3 );
	const Vec t6 = v0 * ( nm3 * v0 + v3 + t4 );
	double tmp = hadd12( t6, v1 * ( nm3 * v1 + v0 + v2 + v3 ), v2 * ( nm3 * v2 + v0 + v1 + v3 ), v3 * ( nm3 * v3 + v0 + t4 ) );
	const double s2 = ( 0.666666666666667 * 0.5 ) * tmp / s0;

	STORE( v0 );
	STORE( v1 );
	STORE( v3 );
	STORE( v2 );
	STORE( t0 );
	STORE( t1 );
	STORE( t2 );
	STORE( t3 );
	STORE( prod );

	result_0[ 0 ] = s1 * ( -v2_x + 3.0 * v0_x + s2 * ( t3_z * t0_y - t3_y * t0_z + t1_z * t2_y - t1_y * t2_z ) - v1_x - v3_x );
	result_0[ 1 ] = s1 * ( 3.0 * v0_y - v3_y - v1_y - v2_y + s2 * ( -prod_y + t0_z * t3_x - t0_x * t3_z ) );
	result_0[ 2 ] = s1 * ( -v1_z - v2_z - v3_z + s2 * ( -t0_y * t3_x - prod_z + t0_x * t3_y ) + 3.0 * v0_z );
}

void floatTetWild::AMIPS_jacobian( const std::array<Scalar, 12>& T, Vector3& result_0 )
{
#if 0
	AMIPS_jacobian_v3( T, result_0 );
#else
	Vector3 resOld, resNew;
	AMIPS_jacobian_v1( T, resOld );
	AMIPS_jacobian_v3( T, resNew );
	Vector3 diff = resNew - resOld;
	result_0 = resOld;
#endif
}