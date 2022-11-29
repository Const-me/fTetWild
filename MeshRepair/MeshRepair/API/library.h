#pragma once
#include "eGlobalFlags.h"
#include "loggerApi.h"

#ifdef _MSC_VER
// On Windows, it's controlled by MeshRepair.def module definition file. There's __declspec(dllexport), but it adds underscore, I don't like that.
#define DLLEXPORT extern "C"
#define COMLIGHTCALL __stdcall
#else
// When compiling the DLL on Linux, assuming users gonna specify `-fvisibility=hidden` compiler flag
#define DLLEXPORT extern "C" __attribute__( ( visibility( "default" ) ) )
#define COMLIGHTCALL
#endif

namespace MeshRepair
{
	struct iMeshRepair;

	DLLEXPORT HRESULT COMLIGHTCALL createMeshRepair( eGlobalFlags globalFlags, const sLoggerSetup* logSetup, iMeshRepair** rdi );
}