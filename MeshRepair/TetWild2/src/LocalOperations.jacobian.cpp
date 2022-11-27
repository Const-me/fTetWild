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

	const double helper_38 = v0_x - v3_x;
	const double helper_3 = v0_y - v3_y;
	const double helper_24 = v0_z - v3_z;

	const double helper_7 = magic1 * v0_x - magic2 * v1_x + magic1 * v3_x;
	const double helper_9 = magic3 * v0_z;
	const double helper_11 = magic3 * v1_z;
	const double helper_13 = magic4 * v2_z;
	const double helper_15 = magic3 * v3_z;
	const double helper_16 = helper_11 - helper_13 + helper_15 + helper_9;
	const double helper_17 = magic1 * v0_z;
	const double helper_18 = magic2 * v1_z;
	const double helper_19 = magic1 * v3_z;
	const double helper_20 = helper_17 - helper_18 + helper_19;
	const double helper_22 = -magic4 * v2_x + magic3 * v0_x + magic3 * v1_x + magic3 * v3_x;
	const double helper_23 = helper_16 * helper_7 - helper_20 * helper_22;
	const double helper_25 = magic3 * v0_y;
	const double helper_27 = magic3 * v1_y;
	const double helper_29 = magic4 * v2_y;
	const double helper_30 = magic3 * v3_y;
	const double helper_31 = helper_25 + helper_27 - helper_29 + helper_30;
	const double helper_32 = helper_31 * helper_7;
	const double helper_33 = magic1 * v0_y;
	const double helper_34 = magic2 * v1_y;
	const double helper_35 = magic1 * v3_y;
	const double helper_36 = helper_33 - helper_34 + helper_35;
	const double helper_37 = helper_22 * helper_36;
	const double helper_39 = helper_23 * helper_3 - helper_24 * ( helper_32 - helper_37 ) - helper_38 * ( helper_16 * helper_36 - helper_20 * helper_31 );
	const double helper_40 = 1.0 / cubicRoot( pow2( helper_39 ) );
	const double helper_41 = magic5 * v1_z - magic5 * v2_z;
	const double helper_42 = magic5 * v1_y - magic5 * v2_y;
	const double helper_43 = 0.5 * v2_x + 0.5 * v1_x;
	const double helper_44 = 0.5 * v1_y + 0.5 * v2_y;
	const double helper_45 = 0.5 * v1_z + 0.5 * v2_z;
	const double helper_46 = 0.666666666666667 *
							 ( v0_y * ( -1.5 * v0_y + 0.5 * v3_y + helper_44 ) + v1_z * ( -1.5 * v1_z + 0.5 * v2_z + 0.5 * v3_z + 0.5 * v0_z ) +
							   v2_z * ( 0.5 * v1_z - 1.5 * v2_z + 0.5 * v3_z + 0.5 * v0_z ) + v3_z * ( -1.5 * v3_z + helper_45 + 0.5 * v0_z ) +
							   v3_y * ( 0.5 * v0_y - 1.5 * v3_y + helper_44 ) + v2_x * ( -1.5 * v2_x + 0.5 * v0_x + 0.5 * v1_x + 0.5 * v3_x ) +
							   v1_y * ( 0.5 * v0_y + 0.5 * v3_y - 1.5 * v1_y + 0.5 * v2_y ) + v2_y * ( 0.5 * v0_y + 0.5 * v3_y + 0.5 * v1_y - 1.5 * v2_y ) +
							   v0_x * ( -1.5 * v0_x + helper_43 + 0.5 * v3_x ) + v1_x * ( 0.5 * v2_x + 0.5 * v0_x - 1.5 * v1_x + 0.5 * v3_x ) +
							   v3_x * ( 0.5 * v0_x + helper_43 - 1.5 * v3_x ) + v0_z * ( 0.5 * v3_z + helper_45 - 1.5 * v0_z ) ) /
							 helper_39;
	const double helper_47 = -magic5 * v2_x + magic5 * v1_x;
	result_0[ 0 ] = -helper_40 * ( v2_x - 3.0 * v0_x +
								   helper_46 * ( helper_41 * ( -v0_y + v3_y ) - helper_42 * ( v3_z - v0_z ) -
												 ( -helper_17 + helper_18 - helper_19 ) * ( -helper_25 - helper_27 + helper_29 - helper_30 ) +
												 ( -helper_33 + helper_34 - helper_35 ) * ( -helper_11 + helper_13 - helper_15 - helper_9 ) ) +
								   v1_x + v3_x );
	result_0[ 1 ] = helper_40 * ( 3.0 * v0_y - v3_y - v1_y - v2_y + helper_46 * ( helper_23 + helper_24 * helper_47 - helper_38 * helper_41 ) );
	result_0[ 2 ] = helper_40 * ( -v1_z - v2_z - v3_z + helper_46 * ( -helper_3 * helper_47 - helper_32 + helper_37 + helper_38 * helper_42 ) + 3.0 * v0_z );
}

void floatTetWild::AMIPS_jacobian( const std::array<Scalar, 12>& T, Vector3& result_0 )
{
#if 0
	AMIPS_jacobian_v2( T, result_0 );
#else
	Vector3 resOld, resNew;
	AMIPS_jacobian_v1( T, resOld );
	AMIPS_jacobian_v2( T, resNew );
	Vector3 diff = resNew - resOld;
	if( diff != Vector3 ::Zero() )
		__debugbreak();
	result_0 = resOld;
#endif
}