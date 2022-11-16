#include "stdafx.h"
#include "TimeMeasure.h"

namespace
{
	inline uint64_t tscNow()
	{
		// QueryPerformanceCounter() on Windows, and clock_gettime() on Linux, are kernel calls.
		// We don't want them for this use case, too expensive.
		// Fortunately, CPUs support RDTSC instruction for that use case.
		// On modern CPUs, that instruction no longer counts cycles, instead it measures the time, based on CPU's base frequency.
		return __rdtsc();
	}
}  // namespace

TimeMeasure::MeasureRaii::MeasureRaii( TimeMeasure& tm )
	: source( tm )
	, started( tscNow() )
{
}

TimeMeasure::MeasureRaii::~MeasureRaii()
{
	uint64_t elapsed = tscNow() - started;
	source.count++;
	source.time += elapsed;
}

namespace
{
	inline uint64_t qpcNow()
	{
		uint64_t res;
		QueryPerformanceCounter( (LARGE_INTEGER*)&res );
		return res;
	}

	inline uint64_t qpcFrequency()
	{
		uint64_t res;
		QueryPerformanceFrequency( (LARGE_INTEGER*)&res );
		return res;
	}

	inline double computeMultiplier( uint64_t constructedTsc, uint64_t constructedQpc, const Logger& log )
	{
		const uint64_t elapsedTsc = tscNow() - constructedTsc;
		const uint64_t elapsedQpc = qpcNow() - constructedQpc;
		const uint64_t qpcFreq = qpcFrequency();
		const double mul = (double)(int64_t)elapsedQpc / ( (double)(int64_t)qpcFreq * (double)(int64_t)elapsedTsc );

		const double GHz = 1E-9 / mul;
		log.logDebug( "Computed CPU base frequency: %g GHz", GHz );
		return mul;
	}

	inline double seconds( uint64_t elapsed, double mul )
	{
		// CPUs support cvtsi2sd instruction to convert int64_t to FP64
		// They don't have an instruction to convert uint64_t to FP64, that's a library function, relatively expensive one
		double e = (double)(int64_t)elapsed;
		return e * mul;
	}

	inline double printSeconds( double s, const char** unit )
	{
		if( s >= 1 || s <= 0 )
		{
			*unit = "seconds";
			return s;
		}
		if( s >= 1E-3 )
		{
			*unit = "milliseconds";
			return s * 1E3;
		}
		*unit = "microseconds";
		return s * 1E6;
	}
}  // namespace

void TimeMeasure::logInfo( const Logger& log, double mulSeconds, const char* what ) const
{
	const size_t calls = count.load();

	double totalTime = seconds( time.load(), mulSeconds );
	double averageTime = 0;
	if( calls > 0 )
		averageTime = totalTime / (double)(int64_t)calls;

	const char *ut, *ua;
	totalTime = printSeconds( totalTime, &ut );
	averageTime = printSeconds( averageTime, &ua );

	log.logInfo( "%s: %g %s average, %zu calls, %g %s total", what, averageTime, ua, calls, totalTime, ut );
}

TimeMeasures::TimeMeasures()
	: constructedTsc( tscNow() )
	, constructedQpc( qpcNow() )
{
}

void TimeMeasures::logInfo( const Logger& log ) const
{
	const double mul = computeMultiplier( constructedTsc, constructedQpc, log );

#define LOG_ENTRY( F ) F.logInfo( log, mul, #F )

	LOG_ENTRY( edgeCollapsingAux );
	LOG_ENTRY( collapseAnEdge );
	LOG_ENTRY( insertTrianglesAux );
	LOG_ENTRY( insertOneTriangle );
	LOG_ENTRY( findCuttingTets );
	LOG_ENTRY( edgeSwapping );
}