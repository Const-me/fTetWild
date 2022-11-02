#include "stdafx.h"
#include "TriangleMesh.h"
#include "setMeshData.h"
#include "BoundingBox.hpp"

// Set vertex buffer of the mesh, upcasting coordinates to FP64
HRESULT TriangleMesh::assignVertices( size_t count, const float* vb )
{
	try
	{
		vertices.resize( count );
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}
	upcastFloats( (double*)vertices.data(), count * 3, vb );
	return S_OK;
}

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