#include "stdafx.h"
#include "MeshBase.h"
#include "../miscUtils.h"

// Set vertex buffer of the mesh, upcasting coordinates to FP64
HRESULT MeshBase::assignVertices( size_t count, const float* vb )
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