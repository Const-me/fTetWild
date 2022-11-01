#pragma once
#include "../../ComLightLib/comLightClient.h"
#include "Parameters.h"

namespace MeshRepair
{
	struct DECLSPEC_NOVTABLE iSourceMesh : public ComLight::IUnknown
	{
		DEFINE_INTERFACE_ID( "{8cb3666b-83ac-4969-b0a2-6f628c7db797}" );
	};

	struct DECLSPEC_NOVTABLE iResultMesh : public ComLight::IUnknown
	{
		DEFINE_INTERFACE_ID( "{b78deaa8-d975-415c-8ec6-e075dfdb8fea}" );

		virtual HRESULT COMLIGHTCALL getSize( uint32_t& vertices, uint32_t& triangles ) = 0;

		// Length written = triangles * 3
		virtual HRESULT COMLIGHTCALL getIndexBuffer( uint32_t* rdi ) = 0;
		// Length written = vertices * 3
		virtual HRESULT COMLIGHTCALL getVertexBufferFP32( float* rdi ) = 0;
		// Length written = vertices * 3
		virtual HRESULT COMLIGHTCALL getVertexBufferFP64( double* rdi ) = 0;
	};

	struct DECLSPEC_NOVTABLE iMeshRepair : public ComLight::IUnknown
	{
		DEFINE_INTERFACE_ID( "{f4ef4dfd-6016-45bb-ac93-2865c9c044e6}" );

		virtual HRESULT COMLIGHTCALL createIndexedMeshFP32(
		  uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib, iSourceMesh** rdi ) = 0;

		virtual HRESULT COMLIGHTCALL repair( iSourceMesh* mesh, const Parameters& parameters, iResultMesh** rdi ) = 0;
	};
}  // namespace MeshRepair