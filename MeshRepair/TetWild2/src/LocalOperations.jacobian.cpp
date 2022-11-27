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
}  // namespace

void floatTetWild::AMIPS_jacobian_v2( const std::array<Scalar, 12>& arr, Vector3& result_0 )
{
	const double helper_4 = arr[ 0 ];
	const double helper_1 = arr[ 1 ];
	const double helper_8 = arr[ 2 ];
	const double helper_5 = arr[ 3 ];
	const double helper_26 = arr[ 4 ];
	const double helper_10 = arr[ 5 ];
	const double helper_21 = arr[ 6 ];
	const double helper_28 = arr[ 7 ];
	const double helper_12 = arr[ 8 ];
	const double helper_6 = arr[ 9 ];
	const double helper_2 = arr[ 10 ];
	const double helper_14 = arr[ 11 ];

	const double helper_3 = helper_1 - helper_2;
	const double helper_7 = 0.577350269189626 * helper_4 - 1.15470053837925 * helper_5 + 0.577350269189626 * helper_6;
	const double helper_9 = 0.408248290463863 * helper_8;
	const double helper_11 = 0.408248290463863 * helper_10;
	const double helper_13 = 1.22474487139159 * helper_12;
	const double helper_15 = 0.408248290463863 * helper_14;
	const double helper_16 = helper_11 - helper_13 + helper_15 + helper_9;
	const double helper_17 = 0.577350269189626 * helper_8;
	const double helper_18 = 1.15470053837925 * helper_10;
	const double helper_19 = 0.577350269189626 * helper_14;
	const double helper_20 = helper_17 - helper_18 + helper_19;
	const double helper_22 = -1.22474487139159 * helper_21 + 0.408248290463863 * helper_4 + 0.408248290463863 * helper_5 + 0.408248290463863 * helper_6;
	const double helper_23 = helper_16 * helper_7 - helper_20 * helper_22;
	const double helper_24 = -helper_14 + helper_8;
	const double helper_25 = 0.408248290463863 * helper_1;
	const double helper_27 = 0.408248290463863 * helper_26;
	const double helper_29 = 1.22474487139159 * helper_28;
	const double helper_30 = 0.408248290463863 * helper_2;
	const double helper_31 = helper_25 + helper_27 - helper_29 + helper_30;
	const double helper_32 = helper_31 * helper_7;
	const double helper_33 = 0.577350269189626 * helper_1;
	const double helper_34 = 1.15470053837925 * helper_26;
	const double helper_35 = 0.577350269189626 * helper_2;
	const double helper_36 = helper_33 - helper_34 + helper_35;
	const double helper_37 = helper_22 * helper_36;
	const double helper_38 = helper_4 - helper_6;
	const double helper_39 = helper_23 * helper_3 - helper_24 * ( helper_32 - helper_37 ) - helper_38 * ( helper_16 * helper_36 - helper_20 * helper_31 );
	const double helper_40 = 1.0 / cubicRoot( pow2( helper_39 ) );
	const double helper_41 = 0.707106781186548 * helper_10 - 0.707106781186548 * helper_12;
	const double helper_42 = 0.707106781186548 * helper_26 - 0.707106781186548 * helper_28;
	const double helper_43 = 0.5 * helper_21 + 0.5 * helper_5;
	const double helper_44 = 0.5 * helper_26 + 0.5 * helper_28;
	const double helper_45 = 0.5 * helper_10 + 0.5 * helper_12;
	const double helper_46 =
	  0.666666666666667 *
	  ( helper_1 * ( -1.5 * helper_1 + 0.5 * helper_2 + helper_44 ) + helper_10 * ( -1.5 * helper_10 + 0.5 * helper_12 + 0.5 * helper_14 + 0.5 * helper_8 ) +
		helper_12 * ( 0.5 * helper_10 - 1.5 * helper_12 + 0.5 * helper_14 + 0.5 * helper_8 ) + helper_14 * ( -1.5 * helper_14 + helper_45 + 0.5 * helper_8 ) +
		helper_2 * ( 0.5 * helper_1 - 1.5 * helper_2 + helper_44 ) + helper_21 * ( -1.5 * helper_21 + 0.5 * helper_4 + 0.5 * helper_5 + 0.5 * helper_6 ) +
		helper_26 * ( 0.5 * helper_1 + 0.5 * helper_2 - 1.5 * helper_26 + 0.5 * helper_28 ) +
		helper_28 * ( 0.5 * helper_1 + 0.5 * helper_2 + 0.5 * helper_26 - 1.5 * helper_28 ) + helper_4 * ( -1.5 * helper_4 + helper_43 + 0.5 * helper_6 ) +
		helper_5 * ( 0.5 * helper_21 + 0.5 * helper_4 - 1.5 * helper_5 + 0.5 * helper_6 ) + helper_6 * ( 0.5 * helper_4 + helper_43 - 1.5 * helper_6 ) +
		helper_8 * ( 0.5 * helper_14 + helper_45 - 1.5 * helper_8 ) ) /
	  helper_39;
	const double helper_47 = -0.707106781186548 * helper_21 + 0.707106781186548 * helper_5;
	result_0[ 0 ] = -helper_40 * ( 1.0 * helper_21 - 3.0 * helper_4 +
								   helper_46 * ( helper_41 * ( -helper_1 + helper_2 ) - helper_42 * ( helper_14 - helper_8 ) -
												 ( -helper_17 + helper_18 - helper_19 ) * ( -helper_25 - helper_27 + helper_29 - helper_30 ) +
												 ( -helper_33 + helper_34 - helper_35 ) * ( -helper_11 + helper_13 - helper_15 - helper_9 ) ) +
								   1.0 * helper_5 + 1.0 * helper_6 );
	result_0[ 1 ] = helper_40 * ( 3.0 * helper_1 - 1.0 * helper_2 - 1.0 * helper_26 - 1.0 * helper_28 +
								  helper_46 * ( helper_23 + helper_24 * helper_47 - helper_38 * helper_41 ) );
	result_0[ 2 ] = helper_40 * ( -1.0 * helper_10 - 1.0 * helper_12 - 1.0 * helper_14 +
								  helper_46 * ( -helper_3 * helper_47 - helper_32 + helper_37 + helper_38 * helper_42 ) + 3.0 * helper_8 );
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