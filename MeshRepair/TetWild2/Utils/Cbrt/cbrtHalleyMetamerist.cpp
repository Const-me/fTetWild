#include "stdafx.h"
#include "cbrt.h"

// https://web.archive.org/web/20131227144655/http://metamerist.com/cbrt/CubeRoot.cpp
namespace
{
	// cube root approximation using bit hack for 64-bit float
	// adapted from Kahan's cbrt
	__forceinline double cbrt_5d( double d )
	{
		const unsigned int B1 = 715094163;
		double t = 0.0;
		unsigned int* pt = (unsigned int*)&t;
		unsigned int* px = (unsigned int*)&d;
		pt[ 1 ] = px[ 1 ] / 3 + B1;
		return t;
	}

	// iterative cube root approximation using Halley's method (double)
	__forceinline double cbrta_halleyd( const double a, const double R )
	{
		const double a3 = a * a * a;
		const double b = a * ( a3 + R + R ) / ( a3 + a3 + R );
		return b;
	}

	// cube root approximation using 3 iterations of Halley's method (double)
	double halley_cbrt3d( double d )
	{
		double a = cbrt_5d( d );
		a = cbrta_halleyd( a, d );
		a = cbrta_halleyd( a, d );
		return cbrta_halleyd( a, d );
	}
}

double LowLevel::cbrtHalleyMetamerist( double val )
{
	return halley_cbrt3d( val );
}