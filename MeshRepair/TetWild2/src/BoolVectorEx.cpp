#include "stdafx.h"
#include "BoolVectorEx.h"
using namespace floatTetWild;

namespace
{
	__forceinline void writeZeros( void* rdi )
	{
		const __m256 zz = _mm256_setzero_ps();
		float* pf = (float*)rdi;
		_mm256_store_ps( pf, zz );
		_mm256_store_ps( pf + 8, zz );
	}
}  // namespace

void BoolVectorEx::initEmpty( size_t len )
{
	if( length > 0 )
	{
		// Clear all the bits previously set
		for( size_t i = 0; i < outer.size(); i++ )
		{
			uint64_t v = outer[ i ];
			if( v == 0 )
				continue;
			outer[ i ] = 0;

			const size_t innerBase = i * 64;
			do
			{
				const size_t bitIndex = _tzcnt_u64( v );
				v = _blsr_u64( v );
				const size_t innerBlock = innerBase + bitIndex;
				writeZeros( &inner[ innerBlock ] );
			} while( 0 != v );
		}

		const size_t prevInner = inner.size();
		const size_t countInner = ( len + 511 ) / 512;
		inner.resize( countInner );
		if( countInner > prevInner )
			memset( &inner[ prevInner ], 0, ( countInner - prevInner ) * 64 );

		const size_t countOuter = ( countInner + 63 ) / 64;
		outer.resize( countOuter, 0 );
		length = len;
	}
	else
	{
		const size_t countInner = ( len + 511 ) / 512;
		inner.resize( countInner );
		memset( inner.data(), 0, countInner * 64 );

		const size_t countOuter = ( countInner + 63 ) / 64;
		outer.resize( countOuter, 0 );

		length = len;
		return;
	}
}

namespace
{
	__forceinline void enqueueBlock( const void* pv, size_t indexBase, std::queue<int>& queue )
	{
		const uint64_t* rsi = (const uint64_t*)pv;
		size_t indexLast = indexBase + 512;
		for( ; indexBase < indexLast; indexBase += 64, rsi++ )
		{
			uint64_t v = *rsi;
			if( v == 0 )
				continue;

			do
			{
				const size_t bitIndex = _tzcnt_u64( v );
				v = _blsr_u64( v );
				const size_t idx = indexBase + bitIndex;
				queue.push( (int)idx );
			} while( 0 != v );
		}
	}
}  // namespace

void BoolVectorEx::enqueueSet( std::queue<int>& queue ) const
{
	for( size_t i = 0; i < outer.size(); i++ )
	{
		uint64_t v = outer[ i ];
		if( v == 0 )
			continue;

		const size_t innerBase = i * 64;
		do
		{
			const size_t bitIndex = _tzcnt_u64( v );
			v = _blsr_u64( v );
			const size_t innerBlock = innerBase + bitIndex;

			enqueueBlock( &inner[ innerBlock ], innerBlock * 512, queue );
		} while( 0 != v );
	}
}