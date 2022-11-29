// Windows native flavor of the COM interfaces implemented by the MeshRepair.dll
// Include this header if you're consuming the library from a Win32 application or DLL
#pragma once
#ifndef _MSC_VER
#error This header is only useful on Windows when compiling with MSVC
#endif

// Include Windows.h, for HRESULT codes
#ifndef _WINDOWS_
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif
// That header defines IUnknown interface
#include <unknwn.h>
// Parameters.h includes <stdint.h>
#include "Parameters.h"

namespace MeshRepair
{
	struct __declspec( novtable, uuid( "8cb3666b-83ac-4969-b0a2-6f628c7db797" ) ) iSourceMesh : public IUnknown
	{
	};

	struct __declspec( novtable, uuid( "b78deaa8-d975-415c-8ec6-e075dfdb8fea" ) ) iResultMesh : public IUnknown
	{
		virtual HRESULT __stdcall getSize( uint32_t& vertices, uint32_t& triangles ) const = 0;

		// Length written = triangles * 3
		virtual HRESULT __stdcall getIndexBuffer( uint32_t* rdi ) const = 0;
		// Length written = vertices * 3
		virtual HRESULT __stdcall getVertexBufferFP32( float* rdi ) const = 0;
		// Length written = vertices * 3
		virtual HRESULT __stdcall getVertexBufferFP64( double* rdi ) const = 0;
	};

	struct __declspec( novtable, uuid( "f4ef4dfd-6016-45bb-ac93-2865c9c044e6" ) ) iMeshRepair : public IUnknown
	{
		virtual HRESULT __stdcall createIndexedMeshFP32(
		  uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib, iSourceMesh** rdi ) = 0;

		virtual HRESULT __stdcall repair( iSourceMesh* mesh, const Parameters& parameters, iResultMesh** rdi ) = 0;
	};
}  // namespace MeshRepair