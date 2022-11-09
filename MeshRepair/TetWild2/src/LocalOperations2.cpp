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
	const double a2 = arr[ 2 ];
	const double a11 = arr[ 11 ];
	const double a0 = arr[ 0 ];
	const double a3 = arr[ 3 ];
	const double a9 = arr[ 9 ];
	const double tmp0 = magic1 * a0 - magic2 * a3 + magic1 * a9;
	const double a1 = arr[ 1 ];
	const double a4 = arr[ 4 ];
	const double a7 = arr[ 7 ];
	const double a10 = arr[ 10 ];
	const double tmp1 = magic3 * a10 + magic3 * a1 + magic3 * a4 - magic4 * a7;
	const double tmp2 = magic1 * a10 + magic1 * a1 - magic2 * a4;
	const double a6 = arr[ 6 ];
	const double tmp3 = magic3 * a0 - magic4 * a6+ magic3 * a3 + magic3 * a9;
	const double a5 = arr[ 5 ];
	const double a8 = arr[ 8 ];
	const double tmp4 = magic3 * a2 + magic3 * a5 + magic3 * a11 - magic4 * a8;
	const double tmp5 = magic1 * a2 + magic1 * a11 - magic2 * a5 ;

	const double tmp6 = 0.5 * ( a6 + a3 );
	const double tmp7 = 0.5 * ( a4 + a7 );
	const double tmp8 = 0.5 * ( a5 + a8 );

	const double tmpDiv = 
		( a2 - a11 ) * ( tmp1 * tmp0 - tmp2 * tmp3 ) -
		( a1 - a10 ) * ( tmp4 * tmp0 - tmp3 * tmp5 ) +
		( a0 - a9  ) * ( tmp2 * tmp4 - tmp1 * tmp5  );

	const double tmpMul = 
		a2 * ( -1.5 * a2 + 0.5 * a11 + tmp8 ) +
		a10 * ( -1.5 * a10 + tmp7 + 0.5 * a1 ) + 
		a6 * ( -1.5 * a6 + 0.5 * a0 + 0.5 * a3 + 0.5 * a9 ) +
		a5 * ( 0.5 * a2 - 1.5 * a5 + 0.5 * a8 + 0.5 * a11 ) + 
		a8 * ( 0.5 * a2 + 0.5 * a5 - 1.5 * a8 + 0.5 * a11 ) +
		a11 * ( 0.5 * a2 - 1.5 * a11 + tmp8 ) + 
		a0 * ( tmp6 - 1.5 * a0 + 0.5 * a9 ) + 
		a3 * ( 0.5 * a6 + 0.5 * a0 - 1.5 * a3 + 0.5 * a9 ) +
		a9 * ( tmp6 + 0.5 * a0 - 1.5 * a9 ) + 
		a1 * ( 0.5 * a10 + tmp7 - 1.5 * a1 ) + 
		a4 * ( 0.5 * a10 + 0.5 * a1 - 1.5 * a4 + 0.5 * a7 ) +
		a7 * ( 0.5 * a10 + 0.5 * a1 + 0.5 * a4 - 1.5 * a7 );
	
	double res = -( tmpMul ) / std::cbrt( tmpDiv * tmpDiv );

	return res;
}