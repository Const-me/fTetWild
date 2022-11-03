#include "stdafx.h"
#include "orient3D.h"
#include "RobustPredicates.h"

namespace GEO2
{
	constexpr int FPG_UNCERTAIN_VALUE = 0;

	static Sign orient_3d_filter( const double* p0, const double* p1, const double* p2, const double* p3 )
	{
		double a11;
		a11 = ( p1[ 0 ] - p0[ 0 ] );
		double a12;
		a12 = ( p1[ 1 ] - p0[ 1 ] );
		double a13;
		a13 = ( p1[ 2 ] - p0[ 2 ] );
		double a21;
		a21 = ( p2[ 0 ] - p0[ 0 ] );
		double a22;
		a22 = ( p2[ 1 ] - p0[ 1 ] );
		double a23;
		a23 = ( p2[ 2 ] - p0[ 2 ] );
		double a31;
		a31 = ( p3[ 0 ] - p0[ 0 ] );
		double a32;
		a32 = ( p3[ 1 ] - p0[ 1 ] );
		double a33;
		a33 = ( p3[ 2 ] - p0[ 2 ] );
		double Delta;
		Delta = ( ( ( a11 * ( ( a22 * a33 ) - ( a23 * a32 ) ) ) - ( a21 * ( ( a12 * a33 ) - ( a13 * a32 ) ) ) ) + ( a31 * ( ( a12 * a23 ) - ( a13 * a22 ) ) ) );
		double eps;
		double max1 = fabs( a11 );
		if( ( max1 < fabs( a21 ) ) )
			max1 = fabs( a21 );
		if( ( max1 < fabs( a31 ) ) )
			max1 = fabs( a31 );
		double max2 = fabs( a12 );
		if( ( max2 < fabs( a13 ) ) )
			max2 = fabs( a13 );
		if( ( max2 < fabs( a22 ) ) )
			max2 = fabs( a22 );
		if( ( max2 < fabs( a23 ) ) )
			max2 = fabs( a23 );
		double max3 = fabs( a22 );
		if( ( max3 < fabs( a23 ) ) )
			max3 = fabs( a23 );
		if( ( max3 < fabs( a32 ) ) )
			max3 = fabs( a32 );
		if( ( max3 < fabs( a33 ) ) )
			max3 = fabs( a33 );
		double lower_bound_1;
		double upper_bound_1;
		lower_bound_1 = max1;
		upper_bound_1 = max1;
		if( ( max2 < lower_bound_1 ) )
			lower_bound_1 = max2;
		else if( ( max2 > upper_bound_1 ) )
			upper_bound_1 = max2;
		if( ( max3 < lower_bound_1 ) )
			lower_bound_1 = max3;
		else if( ( max3 > upper_bound_1 ) )
			upper_bound_1 = max3;

		if( ( lower_bound_1 < 1.63288018496748314939e-98 ) )
			return Sign::ZERO;
		if( ( upper_bound_1 > 5.59936185544450928309e+101 ) )
			return Sign::ZERO;

		eps = ( 5.11071278299732992696e-15 * ( ( max2 * max3 ) * max1 ) );
		if( ( Delta > eps ) )
			return Sign::POSITIVE;
		if( ( Delta < -eps ) )
			return Sign::NEGATIVE;
		return Sign::ZERO;
	}

	Sign orient_3d( const double* p0, const double* p1, const double* p2, const double* p3 )
	{
		double d = RobustPredicates::orient3d( p0, p1, p2, p3 );
		Sign result = orient_3d_filter( p0, p1, p2, p3 );
		__debugbreak();
		return result;
	}
}  // namespace GEO2