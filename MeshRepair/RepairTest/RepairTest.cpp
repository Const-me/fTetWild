#include "stdafx.h"
#include "Utils/IndexedMesh.h"
#include "../MeshRepair/API/library.h"
#include "Utils/Timer.h"
#include "Utils/ConsoleLogSink.h"
#include "Utils/formatMessage.h"

// static const LPCTSTR stlSource = LR"(C:\Temp\2remove\MeshRepair\model.stl)";
static const LPCTSTR stlSource = LR"(C:\Temp\2remove\MeshRepair\dragon.stl)";

static CString resultPath( LPCTSTR source )
{
	LPCTSTR ext = PathFindExtensionW( source );
	if( nullptr == ext || 0 == *ext )
	{
		fprintf( stderr, "Please supply a complete file name, including the .stl extension\n" );
		throw E_INVALIDARG;
	}

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

HRESULT repairMesh( LPCTSTR stl, CString& res )
{
	IndexedMesh mesh;
	CHECK( mesh.loadBinaryStl( stl ) );
	res = resultPath( stl );

	using ComLight::CComPtr;
	using namespace MeshRepair;
	CComPtr<iMeshRepair> repair;
	constexpr eLogLevel consoleLogLevel = eLogLevel::Debug;
	ConsoleLogSink consoleLogSink { consoleLogLevel };
	const eGlobalFlags globalFlags = eGlobalFlags::UseOpenMP;
	CHECK( createMeshRepair( globalFlags, consoleLogSink, &repair ) );

	CComPtr<iSourceMesh> source;
	CHECK( createMesh( repair, mesh, source ) );

	// The defaults are reasonable
	Parameters params;
	// params.flags |= eRepairFlags::SkipSimplify;

	// The default epsilon is 1E-3
	// With 1E-4 it takes about 13 minutes, for the dragon https://en.wikipedia.org/wiki/Stanford_dragon
	// params.epsilon = 1E-4;

	CComPtr<iResultMesh> result;
	{
		Timer tt { "iMeshRepair.repair()" };
		CHECK( repair->repair( source, params, &result ) );
	}

	CHECK( copyMesh( result, mesh ) );
	CHECK( mesh.saveBinaryStl( res ) );
	return S_OK;
}

int wmain( int argc, wchar_t* argv[] )
{
	if( argc != 2 )
	{
		printf( "Usage: RepairTest input.stl\n" );
		return 1;
	}
	HRESULT hr = E_UNEXPECTED;
	CString result;
	try
	{
		hr = repairMesh( argv[ 1 ], result );
	}
	catch( HRESULT h )
	{
		hr = h;
	}
	if( SUCCEEDED( hr ) )
	{
		wprintf( L"Repaired the mesh. The output saved to %s\n", cstr( result ) );
		return 0;
	}

	CString msg = formatMessage( hr );
	fwprintf( stderr, L"Mesh repair failed: %s\n", cstr( msg ) );
	return (int)hr;
}