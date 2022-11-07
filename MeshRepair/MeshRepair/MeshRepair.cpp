#include "stdafx.h"
#include "API/library.h"
#include "SourceMesh.h"
#include <geogram/basic/common.h>
#include "meshRepairMain.h"
// #include "../TetWild2/Utils/RobustPredicates.h"
#include "../GeogramDelaunay/geogram/numerics/predicates.h"

namespace MeshRepair
{
	class MeshRepair : public ComLight::ObjectRoot<iMeshRepair>
	{
		HRESULT COMLIGHTCALL createIndexedMeshFP32(
		  uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib, iSourceMesh** rdi ) noexcept override final;

		HRESULT COMLIGHTCALL repair( iSourceMesh* mesh, const Parameters& parameters, iResultMesh** rdi ) noexcept override final;

	  protected:
		HRESULT FinalConstruct();
		void FinalRelease();
	};

	HRESULT MeshRepair::FinalConstruct()
	{
		// GEO::initialize();
		GEO::PCK::initialize();
		return S_OK;
	}

	void MeshRepair::FinalRelease()
	{
		// GEO::terminate();
		GEO::PCK::terminate();
	}

	DLLEXPORT HRESULT COMLIGHTCALL createMeshRepair( iMeshRepair** rdi )
	{
		// RobustPredicates::exactinit();
		return ComLight::Object<MeshRepair>::create( rdi );
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
			return meshRepairMain( *objMesh, parameters, rdi );
		}
		catch( const std::bad_alloc& )
		{
			return E_OUTOFMEMORY;
		}
	}
}  // namespace MeshRepair