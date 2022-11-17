#include "stdafx.h"
#include "BoolVector.h"
using namespace floatTetWild;

uint64_t BoolVector::countTrue() const
{
	if( !vec.empty() )
	{
		const uint64_t* rsi = vec.data();
		const size_t len = vec.size();
		const uint64_t* const rsiEnd = rsi + len;

		// AMD Zen3 can compute 4x 64-bit POPCNT every clock, and 4x integer additions.
		// Using 4 independent accumulators to saturate throughput instead of stalling on latency
		uint64_t a0 = 0;
		uint64_t a1 = 0;
		uint64_t a2 = 0;
		uint64_t a3 = 0;

		const uint64_t* const rsiEndAligned = rsi + ( len / 4 ) * 4;
		for( ; rsi < rsiEndAligned; rsi += 4 )
		{
			a0 += __popcnt64( rsi[ 0 ] );
			a1 += __popcnt64( rsi[ 1 ] );
			a2 += __popcnt64( rsi[ 2 ] );
			a3 += __popcnt64( rsi[ 3 ] );
		}

		// Reduce to a single accumulator in a0
		a0 += a1;
		a2 += a3;
		a0 += a2;

		// Handle the remainder
		for( ; rsi < rsiEnd; rsi++ )
			a0 += __popcnt64( *rsi );
		return a0;
	}
	return 0;
}

void BoolVector::resize( size_t len, bool fillValue )
{
	if( len == length )
		return;

	if( size() == 0 )
	{
		length = len;
		uint64_t fv = 0;
		if( fillValue )
			fv = ~fv;
		vec.resize( ( len + 63 ) / 64, fv );
		const uint64_t rem = len % 64;
		if( 0 == rem )
			return;
		uint64_t leftover = 64 - rem;
		fv >>= leftover;
		vec.back() = fv;
		return;
	}

	throw std::exception( "BoolVector dynamic resizing is not implemented" );
}