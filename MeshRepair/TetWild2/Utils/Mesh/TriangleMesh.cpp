#include "stdafx.h"
#include "../Geogram2.h"
#include "TriangleMesh.h"
#include "../BoundingBox.hpp"

HRESULT TriangleMesh::assignTriangles( size_t count, const uint32_t* ib )
{
	static_assert( sizeof( vec3i ) == 12, "sizeof" );

	try
	{
		const vec3i* rsi = (const vec3i*)ib;
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
	for( const vec3& vert : vertices )
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

void TriangleMesh::copyDataCpp( std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib ) const
{
	static_assert( sizeof( floatTetWild::Vector3 ) == sizeof( GEO2::vec3 ), "sizeof" );
	static_assert( sizeof( floatTetWild::Vector3i ) == sizeof( GEO2::vec3i ), "sizeof" );

	// Extract the data, store in different types
	const floatTetWild::Vector3* vertBegin = (const floatTetWild::Vector3*)vertexPointer();
	vb.assign( vertBegin, vertBegin + countVertices() );

	const floatTetWild::Vector3i* idxBuffer = (const floatTetWild::Vector3i*)trianglePointer();
	ib.assign( idxBuffer, idxBuffer + countTriangles() );
}

HRESULT TriangleMesh::copyData( std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib ) const
{
	try
	{
		copyDataCpp( vb, ib );
		return S_OK;
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}
}