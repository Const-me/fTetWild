#include "stdafx.h"
#include "API/library.h"
#include "../ComLightLib/comLightServer.h"
#include "SourceMesh.h"
#include <geogram/basic/common.h>
#include "meshRepairMain.h"
#include "../GeogramDelaunay/Predicates/RobustPredicates.h"
// #include "../TetWild2/Utils/cutTableData.h"
#include "../GeogramDelaunay/geogram/numerics/predicates.h"
// #include "../TetWild2/Utils/Geogram2.h"
// #include "../TetWild2/Utils/lowLevel.h"
#include "../TetWild2/parallelThreadsImpl.h"
#include "../TetWild2/Utils/Logger.h"

namespace MeshRepair
{
	class MeshRepair : public ComLight::ObjectRoot<iMeshRepair>
	{
		HRESULT COMLIGHTCALL createIndexedMeshFP32(
		  uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib, iSourceMesh** rdi ) noexcept override final;

		HRESULT COMLIGHTCALL repair( iSourceMesh* mesh, const Parameters& parameters, iResultMesh** rdi ) noexcept override final;

		eGlobalFlags globalFlags;
		sLoggerSetup logger;

	  protected:
		HRESULT FinalConstruct();
		void FinalRelease();

	  public:
		void setFlags( eGlobalFlags flags )
		{
			globalFlags = flags;
		}
		void initLogger( const sLoggerSetup* src );
	};

	HRESULT MeshRepair::FinalConstruct()
	{
		// RobustPredicates::exactinit();
		GEO::PCK::initialize();
		// printCutTableData();
		// validateCutTableData();
		// GEO2::dbgRunSomeTests();
		// floatTetWild::testLowLevelRoutines();
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

	DLLEXPORT HRESULT COMLIGHTCALL createMeshRepair( eGlobalFlags globalFlags, const sLoggerSetup* logSetup, iMeshRepair** rdi )
	{
		ComLight::CComPtr<ComLight::Object<MeshRepair>> result;
		CHECK( ComLight::Object<MeshRepair>::create( result ) );
		result->setFlags( globalFlags );
		result->initLogger( logSetup );
		result.detach( rdi );
		return S_OK;
	}

	HRESULT COMLIGHTCALL MeshRepair::createIndexedMeshFP32(
	  uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib, iSourceMesh** rdi ) noexcept
	{
		if( nullptr == rdi )
			return E_POINTER;
		Logger logger { this->logger };
		logger.logInfo( "iMeshRepair.createIndexedMeshFP32" );

		using namespace ComLight;
		CComPtr<Object<SourceMesh>> result;
		CHECK( Object<SourceMesh>::create( result ) );
		SetThreadsCountRaii raii( globalFlags );
		CHECK( result->createMesh( countVertices, vb, countTriangles, ib ) );
		result.detach( rdi );
		logger.logInfo( "Created an input mesh" );
		return S_OK;
	}

	HRESULT COMLIGHTCALL MeshRepair::repair( iSourceMesh* mesh, const Parameters& parameters, iResultMesh** rdi ) noexcept
	{
		if( nullptr == mesh || nullptr == rdi )
			return E_POINTER;

		SourceMesh* const objMesh = static_cast<SourceMesh*>( mesh );
		try
		{
			return meshRepairMain( *objMesh, globalFlags, parameters, logger, rdi );
		}
		catch( const std::bad_alloc& )
		{
			return E_OUTOFMEMORY;
		}
	}
}  // namespace MeshRepair