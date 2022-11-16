#include "stdafx.h"
#include "VertexConnectedTets.h"
#include "LocalOperations.h"

void floatTetWild::VertexConnectedTets::ensureSorted() const
{
	if( isSorted )
		return;
	std::sort( vec.begin(), vec.end() );
	isSorted = true;
}

void floatTetWild::setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result )
{
	a.ensureSorted();
	b.ensureSorted();
	std::set_intersection( a.vec.begin(), a.vec.end(), b.vec.begin(), b.vec.end(), std::back_inserter( result ) );
}

namespace
{
	static void setIntersectionX3( const int* a1, const int* a2, const int* b1, const int* b2, const int* c1, const int* c2, std::vector<int>& result )
	{
		int a = *a1, b = *b1, c = *c1;
		while( true )
		{
			if( a < b || a < c )
			{
				a1++;
				if( a1 >= a2 )
					return;
				a = *a1;
				continue;
			}
			if( b < a || b < c )
			{
				b1++;
				if( b1 >= b2 )
					return;
				b = *b1;
				continue;
			}
			if( c < a || c < b )
			{
				c1++;
				if( c1 >= c2 )
					return;
				c = *c1;
				continue;
			}
			result.push_back( a );
			a1++;
			b1++;
			c1++;
			if( a1 >= a2 || b1 >= b2 || c1 >= c2 )
				return;
			a = *a1;
			b = *b1;
			c = *c1;
		}
	}
}  // namespace

void floatTetWild::setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result )
{
#if USE_ORIGINAL_VERSION
	set_intersection( a.vec, b.vec, c.vec, result );
#else
	if( a.empty() || b.empty() || c.empty() )
		return;

	a.ensureSorted();
	b.ensureSorted();
	c.ensureSorted();

	const int* p1 = a.vec.data();
	const int* p1End = p1 + a.vec.size();

	const int* p2 = b.vec.data();
	const int* p2End = p2 + b.vec.size();

	const int* p3 = c.vec.data();
	const int* p3End = p3 + c.vec.size();

	setIntersectionX3( p1, p1End, p2, p2End, p3, p3End, result );

#if COMPARE_VERSION
	std::vector<int> orig;
	set_intersection( a.vec, b.vec, c.vec, orig );
	if( result != orig )
		__debugbreak();
#endif
#endif
}