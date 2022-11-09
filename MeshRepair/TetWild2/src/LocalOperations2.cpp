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
	const double a0 = arr[ 0 ];
	const double a1 = arr[ 1 ];
	const double a2 = arr[ 2 ];
	const double a3 = arr[ 3 ];
	const double a4 = arr[ 4 ];
	const double a5 = arr[ 5 ];
	const double a6 = arr[ 6 ];
	const double a7 = arr[ 7 ];
	const double a8 = arr[ 8 ];
	const double a9 = arr[ 9 ];
	const double a10 = arr[ 10 ];
	const double a11 = arr[ 11 ];

	const double helper_6 = 0.577350269189626 * a0 - 1.15470053837925 * a3 + 0.577350269189626 * a9;
	const double helper_11 = 0.408248290463863 * a10 + 0.408248290463863 * a1 + 0.408248290463863 * a4 - 1.22474487139159 * a7;
	const double helper_12 = 0.577350269189626 * a10 + 0.577350269189626 * a1 - 1.15470053837925 * a4;
	const double helper_14 = -1.22474487139159 * a6 + 0.408248290463863 * a0 + 0.408248290463863 * a3 + 0.408248290463863 * a9;

	const double helper_17 = 0.408248290463863 * a2 + 0.408248290463863 * a5 - 1.22474487139159 * a8 + 0.408248290463863 * a11;
	const double helper_18 = 0.577350269189626 * a2 - 1.15470053837925 * a5 + 0.577350269189626 * a11;

	const double helper_19 = 0.5 * a6 + 0.5 * a3;
	const double helper_20 = 0.5 * a4 + 0.5 * a7;
	const double helper_21 = 0.5 * a5 + 0.5 * a8;

	const double tmpDiv = ( a2 - a11 ) * ( helper_11 * helper_6 - helper_12 * helper_14 ) -
					   ( -a10 + a1 ) * ( -helper_14 * helper_18 + helper_17 * helper_6 ) +
					   ( a0 - a9 ) * ( -helper_11 * helper_18 + helper_12 * helper_17 );

	const double res =
	  -( a2 * ( -1.5 * a2 + 0.5 * a11 + helper_21 ) + a10 * ( -1.5 * a10 + helper_20 + 0.5 * a1 ) +
		 a6 * ( -1.5 * a6 + 0.5 * a0 + 0.5 * a3 + 0.5 * a9 ) +
		 a5 * ( 0.5 * a2 - 1.5 * a5 + 0.5 * a8 + 0.5 * a11 ) +
		 a8 * ( 0.5 * a2 + 0.5 * a5 - 1.5 * a8 + 0.5 * a11 ) + a11 * ( 0.5 * a2 - 1.5 * a11 + helper_21 ) +
		 a0 * ( helper_19 - 1.5 * a0 + 0.5 * a9 ) + a3 * ( 0.5 * a6 + 0.5 * a0 - 1.5 * a3 + 0.5 * a9 ) +
		 a9 * ( helper_19 + 0.5 * a0 - 1.5 * a9 ) + a1 * ( 0.5 * a10 + helper_20 - 1.5 * a1 ) +
		 a4 * ( 0.5 * a10 + 0.5 * a1 - 1.5 * a4 + 0.5 * a7 ) +
		 a7 * ( 0.5 * a10 + 0.5 * a1 + 0.5 * a4 - 1.5 * a7 ) ) /
	  std::cbrt( tmpDiv * tmpDiv );

	return res;
}