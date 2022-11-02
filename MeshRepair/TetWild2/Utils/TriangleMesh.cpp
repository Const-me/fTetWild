#include "stdafx.h"
#include "TriangleMesh.h"
#include "setMeshData.h"

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

