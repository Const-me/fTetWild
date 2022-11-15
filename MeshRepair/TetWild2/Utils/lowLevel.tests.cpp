#include "stdafx.h"
#include "lowLevel.h"
#include <atlstr.h>

namespace
{
	static void testCbrt()
	{
		const std::vector<double> testVals = { -4.0, 0.5, 0, -0.0, DBL_MIN, FLT_MIN, 1.0, 8, 15, 378, 1e10 };

		using pfnCbrt = double(*)(double);
		std::vector<double> stdlib, bsdOrig, bsdSse, meta, custom;

		for( double d : testVals )
		{
			stdlib.push_back( std::cbrt( d ) );
			bsdOrig.push_back( LowLevel::cbrtBsdOriginal( d ) );
			bsdSse.push_back( LowLevel::cbrtBsdSse( d ) );
			meta.push_back( LowLevel::cbrtHalleyMetamerist( d ) );
			custom.push_back( LowLevel::cbrtCustom( d ) );
		};

		if( bsdOrig != bsdSse )
			__debugbreak();
		// The rest of them are actually different
		// stdlib has slightly worse precision than BSD, Halley's method fails for very small numbers, Metamerist' version didn't implement subnormals

		__debugbreak();
		CStringA message;
		message.Format( "std %.12g, BSD orig %.12g, BSD SSE %.12g, meta %.12g, custom %.12g\n", stdlib[ 0 ], bsdOrig[ 0 ], bsdSse[ 0 ], meta[ 0 ], custom[ 0 ] );
		OutputDebugStringA( message );
		__debugbreak();
	}
}  // namespace

namespace floatTetWild
{
	void testLowLevelRoutines()
	{
		testCbrt();
	}
}  // namespace floatTetWild