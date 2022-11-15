#include "Cbrt/cbrt.h"

namespace floatTetWild
{
	inline double cbrt( double x )
	{
		return LowLevel::cbrtBsdSse( x );
	}

	void testLowLevelRoutines();
}  // namespace floatTetWild