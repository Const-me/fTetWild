#pragma once

namespace LowLevel
{
	// Copy-pasted from s_cbrt.c, there https://github.com/freebsd/freebsd-src/blob/master/lib/msun/src/s_cbrt.c
	double cbrtBsdOriginal( double x );

	// Ported the above to SSE2. This version appears to be the best of them overall.
	// Sometimes even more precise than standard library, I have used https://www.wolframalpha.com/ as a source of truth
	double cbrtBsdSse( double x );

	// https://web.archive.org/web/20131227144655/http://metamerist.com/cbrt/cbrt.htm
	double cbrtHalleyMetamerist( double x );

	// Incomplete, ain't good
	double cbrtCustom( double val );
}