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

void floatTetWild::AMIPS_jacobian_v2( const std::array<Scalar, 12>& T, Vector3& result_0 )
{
	Scalar helper_0[ 12 ];
	helper_0[ 0 ] = T[ 0 ];
	helper_0[ 1 ] = T[ 1 ];
	helper_0[ 2 ] = T[ 2 ];
	helper_0[ 3 ] = T[ 3 ];
	helper_0[ 4 ] = T[ 4 ];
	helper_0[ 5 ] = T[ 5 ];
	helper_0[ 6 ] = T[ 6 ];
	helper_0[ 7 ] = T[ 7 ];
	helper_0[ 8 ] = T[ 8 ];
	helper_0[ 9 ] = T[ 9 ];
	helper_0[ 10 ] = T[ 10 ];
	helper_0[ 11 ] = T[ 11 ];
	Scalar helper_1 = helper_0[ 1 ];
	Scalar helper_2 = helper_0[ 10 ];
	Scalar helper_3 = helper_1 - helper_2;
	Scalar helper_4 = helper_0[ 0 ];
	Scalar helper_5 = helper_0[ 3 ];
	Scalar helper_6 = helper_0[ 9 ];
	Scalar helper_7 = 0.577350269189626 * helper_4 - 1.15470053837925 * helper_5 + 0.577350269189626 * helper_6;
	Scalar helper_8 = helper_0[ 2 ];
	Scalar helper_9 = 0.408248290463863 * helper_8;
	Scalar helper_10 = helper_0[ 5 ];
	Scalar helper_11 = 0.408248290463863 * helper_10;
	Scalar helper_12 = helper_0[ 8 ];
	Scalar helper_13 = 1.22474487139159 * helper_12;
	Scalar helper_14 = helper_0[ 11 ];
	Scalar helper_15 = 0.408248290463863 * helper_14;
	Scalar helper_16 = helper_11 - helper_13 + helper_15 + helper_9;
	Scalar helper_17 = 0.577350269189626 * helper_8;
	Scalar helper_18 = 1.15470053837925 * helper_10;
	Scalar helper_19 = 0.577350269189626 * helper_14;
	Scalar helper_20 = helper_17 - helper_18 + helper_19;
	Scalar helper_21 = helper_0[ 6 ];
	Scalar helper_22 = -1.22474487139159 * helper_21 + 0.408248290463863 * helper_4 + 0.408248290463863 * helper_5 + 0.408248290463863 * helper_6;
	Scalar helper_23 = helper_16 * helper_7 - helper_20 * helper_22;
	Scalar helper_24 = -helper_14 + helper_8;
	Scalar helper_25 = 0.408248290463863 * helper_1;
	Scalar helper_26 = helper_0[ 4 ];
	Scalar helper_27 = 0.408248290463863 * helper_26;
	Scalar helper_28 = helper_0[ 7 ];
	Scalar helper_29 = 1.22474487139159 * helper_28;
	Scalar helper_30 = 0.408248290463863 * helper_2;
	Scalar helper_31 = helper_25 + helper_27 - helper_29 + helper_30;
	Scalar helper_32 = helper_31 * helper_7;
	Scalar helper_33 = 0.577350269189626 * helper_1;
	Scalar helper_34 = 1.15470053837925 * helper_26;
	Scalar helper_35 = 0.577350269189626 * helper_2;
	Scalar helper_36 = helper_33 - helper_34 + helper_35;
	Scalar helper_37 = helper_22 * helper_36;
	Scalar helper_38 = helper_4 - helper_6;
	Scalar helper_39 = helper_23 * helper_3 - helper_24 * ( helper_32 - helper_37 ) - helper_38 * ( helper_16 * helper_36 - helper_20 * helper_31 );
	Scalar helper_40 = 1.0 / cubicRoot( pow2( helper_39 ) );
	Scalar helper_41 = 0.707106781186548 * helper_10 - 0.707106781186548 * helper_12;
	Scalar helper_42 = 0.707106781186548 * helper_26 - 0.707106781186548 * helper_28;
	Scalar helper_43 = 0.5 * helper_21 + 0.5 * helper_5;
	Scalar helper_44 = 0.5 * helper_26 + 0.5 * helper_28;
	Scalar helper_45 = 0.5 * helper_10 + 0.5 * helper_12;
	Scalar helper_46 =
	  0.666666666666667 *
	  ( helper_1 * ( -1.5 * helper_1 + 0.5 * helper_2 + helper_44 ) + helper_10 * ( -1.5 * helper_10 + 0.5 * helper_12 + 0.5 * helper_14 + 0.5 * helper_8 ) +
		helper_12 * ( 0.5 * helper_10 - 1.5 * helper_12 + 0.5 * helper_14 + 0.5 * helper_8 ) + helper_14 * ( -1.5 * helper_14 + helper_45 + 0.5 * helper_8 ) +
		helper_2 * ( 0.5 * helper_1 - 1.5 * helper_2 + helper_44 ) + helper_21 * ( -1.5 * helper_21 + 0.5 * helper_4 + 0.5 * helper_5 + 0.5 * helper_6 ) +
		helper_26 * ( 0.5 * helper_1 + 0.5 * helper_2 - 1.5 * helper_26 + 0.5 * helper_28 ) +
		helper_28 * ( 0.5 * helper_1 + 0.5 * helper_2 + 0.5 * helper_26 - 1.5 * helper_28 ) + helper_4 * ( -1.5 * helper_4 + helper_43 + 0.5 * helper_6 ) +
		helper_5 * ( 0.5 * helper_21 + 0.5 * helper_4 - 1.5 * helper_5 + 0.5 * helper_6 ) + helper_6 * ( 0.5 * helper_4 + helper_43 - 1.5 * helper_6 ) +
		helper_8 * ( 0.5 * helper_14 + helper_45 - 1.5 * helper_8 ) ) /
	  helper_39;
	Scalar helper_47 = -0.707106781186548 * helper_21 + 0.707106781186548 * helper_5;
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