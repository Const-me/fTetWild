#include "stdafx.h"
#include "API/library.h"
#include "SourceMesh.h"
#include <geogram/basic/common.h>

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
		GEO::initialize();
		return S_OK;
	}

	void MeshRepair::FinalRelease() 
	{
		GEO::terminate();
	}

	DLLEXPORT HRESULT COMLIGHTCALL createMeshRepair( iMeshRepair** rdi )
	{
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
		return E_NOTIMPL;
	}
}  // namespace MeshRepair