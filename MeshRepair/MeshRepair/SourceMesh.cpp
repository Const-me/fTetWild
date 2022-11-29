#include "stdafx.h"
#include "SourceMesh.h"
using namespace MeshRepair;

HRESULT SourceMesh::createMesh( uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib )
{
	if( countVertices == 0 || countTriangles == 0 )
		return E_INVALIDARG;
	if( nullptr == vb || nullptr == ib )
		return E_POINTER;

	CHECK( mesh.assignVertices( countVertices, vb ) );
	CHECK( mesh.assignTriangles( countTriangles, ib ) );

	// See MeshIO::load_mesh function in the original code
	mesh.reorderMorton();

	return S_OK;
}

void SourceMesh::makeBuffers( std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib, std::vector<int>* tags ) const
{
	mesh.copyDataCpp( vb, ib );

	if( nullptr != tags ) 
	{
		tags->clear();
		tags->resize( ib.size(), 0 );
	}
}