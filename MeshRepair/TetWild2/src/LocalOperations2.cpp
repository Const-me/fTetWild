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

	const double helper_6 = magic1 * a0 - magic2 * a3 + magic1 * a9;
	const double helper_11 = magic3 * a10 + magic3 * a1 + magic3 * a4 - magic4 * a7;
	const double helper_12 = magic1 * a10 + magic1 * a1 - magic2 * a4;
	const double helper_14 = -magic4 * a6 + magic3 * a0 + magic3 * a3 + magic3 * a9;

	const double helper_17 = magic3 * a2 + magic3 * a5 - magic4 * a8 + magic3 * a11;
	const double helper_18 = magic1 * a2 - magic2 * a5 + magic1 * a11;

	const double helper_19 = a6 + a3;
	const double helper_20 = a4 + a7;
	const double helper_21 = a5 + a8;

	const double tmpDiv = ( a2 - a11 ) * ( helper_11 * helper_6 - helper_12 * helper_14 ) -
					   ( -a10 + a1 ) * ( -helper_14 * helper_18 + helper_17 * helper_6 ) +
					   ( a0 - a9 ) * ( -helper_11 * helper_18 + helper_12 * helper_17 );

	// This optimization delivers identical results because scaling numbers by powers of 2 (like 0.5) is a lossless operation
	// Applies an offset to exponent, but preserves all these precious mantissa bits
	// That's why it was OK to remove the multiplications by 0.5 in the inner expression
	double mul = a2 * ( -3 * a2 + a11 + helper_21 ) + 
		a10 * ( -3 * a10 + helper_20 + a1 ) +
		a6 * ( -3 * a6 + a0 + a3 + a9 ) +
		a5 * ( a2 - 3 * a5 + a8 + a11 ) +
		a8 * ( a2 + a5 - 3 * a8 + a11 ) + 
		a11 * ( a2 - 3 * a11 + helper_21 ) +
		a0 * ( helper_19 - 3 * a0 + a9 ) + 
		a3 * ( a6 + a0 - 3 * a3 + a9 ) +
		a9 * ( helper_19 + a0 - 3 * a9 ) + 
		a1 * ( a10 + helper_20 - 3 * a1 ) +
		a4 * ( a10 + a1 - 3 * a4 + a7 ) +
		a7 * ( a10 + a1 + a4 - 3 * a7 );

	mul *= -0.5;
	const double res = mul / std::cbrt( tmpDiv * tmpDiv );
	return res;
}

double floatTetWild::AMIPS_energy_aux_v3( const std::array<double, 12>& arr )
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

	const double helper_6 = magic1 * a0 - magic2 * a3 + magic1 * a9;
	const double helper_11 = magic3 * a10 + magic3 * a1 + magic3 * a4 - magic4 * a7;
	const double helper_12 = magic1 * a10 + magic1 * a1 - magic2 * a4;
	const double helper_14 = -magic4 * a6 + magic3 * a0 + magic3 * a3 + magic3 * a9;

	const double helper_17 = magic3 * a2 + magic3 * a5 - magic4 * a8 + magic3 * a11;
	const double helper_18 = magic1 * a2 - magic2 * a5 + magic1 * a11;

	const double tmpDiv = ( a2 - a11 ) * ( helper_11 * helper_6 - helper_12 * helper_14 ) -
					   ( -a10 + a1 ) * ( -helper_14 * helper_18 + helper_17 * helper_6 ) +
					   ( a0 - a9 ) * ( -helper_11 * helper_18 + helper_12 * helper_17 );

	const double sum0 = a0 + a3 + a6 + a9;
	const double sum1 = a1 + a4 + a7 + a10;
	const double sum2 = a2 + a5 + a8 + a11;

	double mul = 
		a0 * ( sum0 - 4 * a0 ) + 
		a1 * ( sum1 - 4 * a1 ) +
		a2 * ( sum2 - 4 * a2 ) + 
		a3 * ( sum0 - 4 * a3 ) +
		a4 * ( sum1 - 4 * a4 ) +
		a5 * ( sum2 - 4 * a5 ) +
		a6 * ( sum0 - 4 * a6 ) +
		a7 * ( sum1 - 4 * a7 ) +
		a8 * ( sum2 - 4 * a8 ) + 
		a9 * ( sum0 - 4 * a9 ) + 
		a10 * ( sum1 - 4 * a10 ) +
		a11 * ( sum2 - 4 * a11 );

	mul *= -0.5;
	const double res = mul / std::cbrt( tmpDiv * tmpDiv );
	return res;
}