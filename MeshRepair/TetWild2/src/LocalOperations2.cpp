#include "stdafx.h"
#include "LocalOperations2.h"
// clang-format off
namespace
{
	constexpr double magic1 = 0.577350269189626;	// sqrt(1/3)
	constexpr double magic2 = 1.15470053837925;		// sqrt(4/3)
	constexpr double magic3 = 0.408248290463863;	// sqrt(1/6)
	constexpr double magic4 = 1.22474487139159;		// sqrt(3/2)
}

double floatTetWild::AMIPS_energy_aux_v2( const std::array<double, 12>& T )
{
	double helper_0[ 12 ];
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

	double helper_1 = helper_0[ 2 ];
	double helper_2 = helper_0[ 11 ];
	double helper_3 = helper_0[ 0 ];
	double helper_4 = helper_0[ 3 ];
	double helper_5 = helper_0[ 9 ];
	double helper_6 = 0.577350269189626 * helper_3 - 1.15470053837925 * helper_4 + 0.577350269189626 * helper_5;
	double helper_7 = helper_0[ 1 ];
	double helper_8 = helper_0[ 4 ];
	double helper_9 = helper_0[ 7 ];
	double helper_10 = helper_0[ 10 ];
	double helper_11 = 0.408248290463863 * helper_10 + 0.408248290463863 * helper_7 + 0.408248290463863 * helper_8 - 1.22474487139159 * helper_9;
	double helper_12 = 0.577350269189626 * helper_10 + 0.577350269189626 * helper_7 - 1.15470053837925 * helper_8;
	double helper_13 = helper_0[ 6 ];
	double helper_14 = -1.22474487139159 * helper_13 + 0.408248290463863 * helper_3 + 0.408248290463863 * helper_4 + 0.408248290463863 * helper_5;
	double helper_15 = helper_0[ 5 ];
	double helper_16 = helper_0[ 8 ];
	double helper_17 = 0.408248290463863 * helper_1 + 0.408248290463863 * helper_15 - 1.22474487139159 * helper_16 + 0.408248290463863 * helper_2;
	double helper_18 = 0.577350269189626 * helper_1 - 1.15470053837925 * helper_15 + 0.577350269189626 * helper_2;
	double helper_19 = 0.5 * helper_13 + 0.5 * helper_4;
	double helper_20 = 0.5 * helper_8 + 0.5 * helper_9;
	double helper_21 = 0.5 * helper_15 + 0.5 * helper_16;
	double helper_22 = ( helper_1 - helper_2 ) * ( helper_11 * helper_6 - helper_12 * helper_14 ) -
					   ( -helper_10 + helper_7 ) * ( -helper_14 * helper_18 + helper_17 * helper_6 ) +
					   ( helper_3 - helper_5 ) * ( -helper_11 * helper_18 + helper_12 * helper_17 );

	double res =
	  -( helper_1 * ( -1.5 * helper_1 + 0.5 * helper_2 + helper_21 ) + helper_10 * ( -1.5 * helper_10 + helper_20 + 0.5 * helper_7 ) +
		 helper_13 * ( -1.5 * helper_13 + 0.5 * helper_3 + 0.5 * helper_4 + 0.5 * helper_5 ) +
		 helper_15 * ( 0.5 * helper_1 - 1.5 * helper_15 + 0.5 * helper_16 + 0.5 * helper_2 ) +
		 helper_16 * ( 0.5 * helper_1 + 0.5 * helper_15 - 1.5 * helper_16 + 0.5 * helper_2 ) + helper_2 * ( 0.5 * helper_1 - 1.5 * helper_2 + helper_21 ) +
		 helper_3 * ( helper_19 - 1.5 * helper_3 + 0.5 * helper_5 ) + helper_4 * ( 0.5 * helper_13 + 0.5 * helper_3 - 1.5 * helper_4 + 0.5 * helper_5 ) +
		 helper_5 * ( helper_19 + 0.5 * helper_3 - 1.5 * helper_5 ) + helper_7 * ( 0.5 * helper_10 + helper_20 - 1.5 * helper_7 ) +
		 helper_8 * ( 0.5 * helper_10 + 0.5 * helper_7 - 1.5 * helper_8 + 0.5 * helper_9 ) +
		 helper_9 * ( 0.5 * helper_10 + 0.5 * helper_7 + 0.5 * helper_8 - 1.5 * helper_9 ) ) /
	  std::cbrt( helper_22 * helper_22 );

	return res;
}