#include "stdafx.h"
#include "SourceMesh.h"
#include "../TetWild2/Utils/setMeshData.h"
#include <geogram/mesh/mesh_reorder.h>
using namespace MeshRepair;

HRESULT SourceMesh::createMesh( uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib )
{
	if( countVertices == 0 || countTriangles == 0 )
		return E_INVALIDARG;
	if( nullptr == vb || nullptr == ib )
		return E_POINTER;

	CHECK( assignMeshVertices( mesh, countVertices, vb ) );
	CHECK( assignMeshTriangles( mesh, countTriangles, ib ) );

	// See MeshIO::load_mesh function in the original code
	GEO::mesh_reorder( mesh, GEO::MESH_ORDER_MORTON );

	// Extract the data, store as different types
	CHECK( copyMeshData( mesh, input_vertices, input_faces ) );
	try
	{
		input_tags.clear();
		input_tags.resize( input_faces.size(), 0 );
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}
	return S_OK;
}