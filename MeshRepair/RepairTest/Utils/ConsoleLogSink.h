#pragma once
#include "../../MeshRepair/API/loggerApi.h"
#include <mutex>

class ConsoleLogSink
{
	MeshRepair::sLoggerSetup setup;
	std::mutex mutex;
	inline void logMessage( MeshRepair::eLogLevel lvl, const char* message );
	static void __cdecl logMessageStatic( void* context, MeshRepair::eLogLevel lvl, const char* message );

	class SetupColors
	{
	  public:
		SetupColors( MeshRepair::eLogLevel level );
		~SetupColors();
	};

  public:
	ConsoleLogSink( MeshRepair::eLogLevel level = MeshRepair::eLogLevel::Info );

	operator const MeshRepair::sLoggerSetup*() const
	{
		return &setup;
	}
};