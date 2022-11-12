#include "stdafx.h"
#include "VertexConnectedTets.h"
#include "LocalOperations.h"

#define USE_ORIGINAL_VERSION 0

namespace
{
	struct alignas( 64 ) IntersectionBuffers
	{
		std::vector<int> a, b, c;
	};
	thread_local IntersectionBuffers buffers;
}  // namespace

void floatTetWild::setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result )
{
#if USE_ORIGINAL_VERSION
	set_intersection( a.vec, b.vec, result );
#else
	IntersectionBuffers& ib = buffers;
	ib.a = a.vec;
	ib.b = b.vec;
	std::sort( ib.a.begin(), ib.a.end() );
	std::sort( ib.b.begin(), ib.b.end() );
	std::set_intersection( ib.a.begin(), ib.a.end(), ib.b.begin(), ib.b.end(), std::back_inserter( result ) );
#endif
}

void floatTetWild::setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result )
{
#if USE_ORIGINAL_VERSION
	set_intersection( a.vec, b.vec, c.vec, result );
#else
	IntersectionBuffers& ib = buffers;
	ib.a = a.vec;
	ib.b = b.vec;
	ib.c = c.vec;
	std::sort( ib.a.begin(), ib.a.end() );
	std::sort( ib.b.begin(), ib.b.end() );
	std::sort( ib.c.begin(), ib.c.end() );

	std::set_intersection( ib.a.begin(), ib.a.end(), ib.b.begin(), ib.b.end(), std::back_inserter( result ) );
	auto it = std::set_intersection( result.begin(), result.end(), ib.c.begin(), ib.c.end(), result.begin() );
	result.resize( it - result.begin() );
#endif
}