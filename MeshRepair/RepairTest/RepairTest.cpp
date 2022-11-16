#include "stdafx.h"
#include "Utils/IndexedMesh.h"
#include "../MeshRepair/API/library.h"
#include "Utils/Timer.h"
#include "Utils/ConsoleLogSink.h"
#include <atlstr.h>

// static const LPCTSTR stlSource = LR"(C:\Temp\2remove\MeshRepair\model.stl)";
static const LPCTSTR stlSource = LR"(C:\Temp\2remove\MeshRepair\dragon.stl)";

static CString resultPath( LPCTSTR source )
{
	LPCTSTR ext = PathFindExtensionW( source );
	if( nullptr == ext || 0 == *ext )
		__debugbreak();
	CString res { source, (int)( ext - source ) };
	res += L"-result.stl";
	return res;
}

namespace MeshRepair
{
	HRESULT createMesh( iMeshRepair* mr, const IndexedMesh& rsi, ComLight::CComPtr<iSourceMesh>& rdi )
	{
		if( rsi.triangles.empty() || rsi.vertices.empty() )
			return OLE_E_BLANK;
		if( rsi.triangles.size() > UINT_MAX || rsi.vertices.size() > UINT_MAX )
			return DISP_E_OVERFLOW;
		rdi.release();
		return mr->createIndexedMeshFP32(
		  (uint32_t)rsi.vertices.size(), (const float*)rsi.vertices.data(), (uint32_t)rsi.triangles.size(), (const uint32_t*)rsi.triangles.data(), &rdi );
	}

	HRESULT copyMesh( iResultMesh* rsi, IndexedMesh& rdi )
	{
		uint32_t verts, tris;
		CHECK( rsi->getSize( verts, tris ) );

		try
		{
#ifdef _DEBUG
			rdi.vertices.clear();
			rdi.triangles.clear();
			DirectX::XMFLOAT3 cv { NAN, NAN, NAN };
			rdi.vertices.resize( verts, cv );
			std::array<uint32_t, 3> cf { 0xCCCCCCCCu, 0xCCCCCCCCu, 0xCCCCCCCCu };
			rdi.triangles.resize( tris, cf );
#else
			rdi.vertices.resize( verts );
			rdi.triangles.resize( tris );
#endif
		}
		catch( const std::bad_alloc& )
		{
			return E_OUTOFMEMORY;
		}

		CHECK( rsi->getVertexBufferFP32( (float*)rdi.vertices.data() ) );
		CHECK( rsi->getIndexBuffer( (uint32_t*)rdi.triangles.data() ) );
		return S_OK;
	}
}  // namespace MeshRepair

HRESULT testStlIO()
{
	IndexedMesh mesh;
	CHECK( mesh.loadBinaryStl( stlSource ) );
	CString res = resultPath( stlSource );
	CHECK( mesh.saveBinaryStl( res ) );
	return S_OK;
}

HRESULT testRepair()
{
	IndexedMesh mesh;
	CHECK( mesh.loadBinaryStl( stlSource ) );
	CString res = resultPath( stlSource );

	using ComLight::CComPtr;
	using namespace MeshRepair;
	CComPtr<iMeshRepair> repair;
	ConsoleLogSink consoleLogSink;
	CHECK( createMeshRepair( consoleLogSink, &repair ) );

	CComPtr<iSourceMesh> source;
	CHECK( createMesh( repair, mesh, source ) );

	Parameters params;
	params.flags |= eRepairFlags::UseOpenMP;

	CComPtr<iResultMesh> result;
	{
		Timer tt { "iMeshRepair.repair()" };
		CHECK( repair->repair( source, params, &result ) );
	}

	CHECK( copyMesh( result, mesh ) );
	CHECK( mesh.saveBinaryStl( res ) );
	return S_OK;
}

int main()
{
	// testStlIO();
	testRepair();
	printf( "Hello World!\n" );
	return 0;
}