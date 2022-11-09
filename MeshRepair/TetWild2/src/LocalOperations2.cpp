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

double floatTetWild::AMIPS_energy_aux_v2( const std::array<double, 12>& arr )
{
	const double helper_1 = arr[ 2 ];
	const double helper_2 = arr[ 11 ];
	const double helper_3 = arr[ 0 ];
	const double helper_4 = arr[ 3 ];
	const double helper_5 = arr[ 9 ];
	const double helper_6 = 0.577350269189626 * helper_3 - 1.15470053837925 * helper_4 + 0.577350269189626 * helper_5;
	const double helper_7 = arr[ 1 ];
	const double helper_8 = arr[ 4 ];
	const double helper_9 = arr[ 7 ];
	const double helper_10 = arr[ 10 ];
	const double helper_11 = 0.408248290463863 * helper_10 + 0.408248290463863 * helper_7 + 0.408248290463863 * helper_8 - 1.22474487139159 * helper_9;
	const double helper_12 = 0.577350269189626 * helper_10 + 0.577350269189626 * helper_7 - 1.15470053837925 * helper_8;
	const double helper_13 = arr[ 6 ];
	const double helper_14 = -1.22474487139159 * helper_13 + 0.408248290463863 * helper_3 + 0.408248290463863 * helper_4 + 0.408248290463863 * helper_5;
	const double helper_15 = arr[ 5 ];
	const double helper_16 = arr[ 8 ];

	const double helper_17 = 0.408248290463863 * helper_1 + 0.408248290463863 * helper_15 - 1.22474487139159 * helper_16 + 0.408248290463863 * helper_2;
	const double helper_18 = 0.577350269189626 * helper_1 - 1.15470053837925 * helper_15 + 0.577350269189626 * helper_2;
	const double helper_19 = 0.5 * helper_13 + 0.5 * helper_4;
	const double helper_20 = 0.5 * helper_8 + 0.5 * helper_9;
	const double helper_21 = 0.5 * helper_15 + 0.5 * helper_16;
	const double helper_22 = ( helper_1 - helper_2 ) * ( helper_11 * helper_6 - helper_12 * helper_14 ) -
					   ( -helper_10 + helper_7 ) * ( -helper_14 * helper_18 + helper_17 * helper_6 ) +
					   ( helper_3 - helper_5 ) * ( -helper_11 * helper_18 + helper_12 * helper_17 );

	const double res =
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