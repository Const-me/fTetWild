#include "stdafx.h"
#include "TriangleMesh.h"
#include "../BoundingBox.hpp"

HRESULT TriangleMesh::assignTriangles( size_t count, const uint32_t* ib )
{
	static_assert( sizeof( GEO::vec3i ) == 12, "sizeof" );

	try
	{
		const GEO::vec3i* rsi = (const GEO::vec3i*)ib;
		triangles.assign( rsi, rsi + count );
		return S_OK;
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}
}

double TriangleMesh::boxDiagonal() const
{
	BoundingBox box;
	for( const GEO::vec3& vert : vertices )
	{
		const double* rsi = (const double*)&vert;
		box.extend( rsi );
	}
	return box.diagonal();
}

void TriangleMesh::clearMesh( bool keepMemory )
{
	vertices.clear();
	triangles.clear();
	if( keepMemory )
		return;

	vertices.shrink_to_fit();
	triangles.shrink_to_fit();
}

HRESULT TriangleMesh::copyData( std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib ) const
{
	static_assert( sizeof( floatTetWild::Vector3 ) == sizeof( GEO::vec3 ), "sizeof" );
	static_assert( sizeof( floatTetWild::Vector3i ) == sizeof( GEO::vec3i ), "sizeof" );

	// Extract the data, store in different types
	try
	{
		const floatTetWild::Vector3* vertBegin = (const floatTetWild::Vector3*)vertexPointer();
		vb.assign( vertBegin, vertBegin + countVertices() );

		const floatTetWild::Vector3i* idxBuffer = (const floatTetWild::Vector3i*)trianglePointer();
		ib.assign( idxBuffer, idxBuffer + countTriangles() );
		return S_OK;
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}
}