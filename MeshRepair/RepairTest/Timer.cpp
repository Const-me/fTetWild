#include "stdafx.h"
#include "Timer.h"

namespace
{
	inline uint64_t counter()
	{
		uint64_t res;
		QueryPerformanceCounter( (LARGE_INTEGER*)&res );
		return res;
	}

	static double multiplierSeconds()
	{
		uint64_t freq;
		QueryPerformanceFrequency( (LARGE_INTEGER*)&freq );
		return 1.0 / (double)(int64_t)( freq );
	}

	static const double mulSeconds = multiplierSeconds();
}  // namespace

Timer::Timer( const char* what )
	: started( counter() )
	, message( what )
{
}

Timer::~Timer()
{
	const uint64_t elapsed = counter() - started;
	double seconds = (double)(int64_t)elapsed;
	seconds *= mulSeconds;
	printf( "%s: %g seconds\n", message, seconds );
}