#include "stdafx.h"
#include "API/library.h"
#include "../ComLightLib/comLightServer.h"
#include "SourceMesh.h"
#include <geogram/basic/common.h>
#include "meshRepairMain.h"
// #include "../TetWild2/Utils/RobustPredicates.h"
// #include "../TetWild2/Utils/cutTableData.h"
#include "../GeogramDelaunay/geogram/numerics/predicates.h"

namespace MeshRepair
{
	class MeshRepair : public ComLight::ObjectRoot<iMeshRepair>
	{
		HRESULT COMLIGHTCALL createIndexedMeshFP32(
		  uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib, iSourceMesh** rdi ) noexcept override final;

		HRESULT COMLIGHTCALL repair( iSourceMesh* mesh, const Parameters& parameters, iResultMesh** rdi ) noexcept override final;

		sLoggerSetup logger;

	  protected:
		HRESULT FinalConstruct();
		void FinalRelease();

	  public:
		void initLogger( const sLoggerSetup* src );
	};

	HRESULT MeshRepair::FinalConstruct()
	{
		GEO::PCK::initialize();
		// printCutTableData();
		// validateCutTableData();
		return S_OK;
	}

	void MeshRepair::FinalRelease()
	{
		GEO::PCK::terminate();
	}

	void MeshRepair::initLogger( const sLoggerSetup* src )
	{
		if( nullptr != src )
			logger = *src;
		else
		{
			logger.sink = nullptr;
			logger.context = nullptr;
			logger.level = eLogLevel::Error;
		}
	}

	DLLEXPORT HRESULT COMLIGHTCALL createMeshRepair( const sLoggerSetup* logSetup, iMeshRepair** rdi )
	{
		ComLight::CComPtr<ComLight::Object<MeshRepair>> result;
		CHECK( ComLight::Object<MeshRepair>::create( result ) );
		result->initLogger( logSetup );
		result.detach( rdi );
		return S_OK;
	}

	HRESULT COMLIGHTCALL MeshRepair::createIndexedMeshFP32(
	  uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib, iSourceMesh** rdi ) noexcept
	{
		if( nullptr == rdi )
			return E_POINTER;

		using namespace ComLight;
		CComPtr<Object<SourceMesh>> result;
		CHECK( Object<SourceMesh>::create( result ) );
		CHECK( result->createMesh( countVertices, vb, countTriangles, ib ) );
		result.detach( rdi );
		return S_OK;
	}

	HRESULT COMLIGHTCALL MeshRepair::repair( iSourceMesh* mesh, const Parameters& parameters, iResultMesh** rdi ) noexcept
	{
		if( nullptr == mesh || nullptr == rdi )
			return E_POINTER;

		SourceMesh* const objMesh = static_cast<SourceMesh*>( mesh );
		try
		{
			return meshRepairMain( *objMesh, parameters, logger, rdi );
		}
		catch( const std::bad_alloc& )
		{
			return E_OUTOFMEMORY;
		}
	}
}  // namespace MeshRepair