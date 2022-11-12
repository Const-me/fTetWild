#include "stdafx.h"
#include "VertexConnectedTets.h"
#include "LocalOperations.h"

void floatTetWild::setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result )
{
	set_intersection( a.vec, b.vec, result );
}

void floatTetWild::setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result )
{
	set_intersection( a.vec, b.vec, c.vec, result );
}