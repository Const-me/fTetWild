#pragma once
#include <stdint.h>

namespace MeshRepair
{
	// Log level for messages
	enum struct eLogLevel : uint8_t
	{
		Error = 0,
		Warning = 1,
		Info = 2,
		Debug = 3
	};

	// C function pointer to receive log messages from the library. The messages are encoded in UTF-8.
	using pfnLoggerSink = void( __cdecl* )( void* context, eLogLevel lvl, const char* message );

	// A sink to receive log messages produced by MeshRepair.dll
	struct sLoggerSetup
	{
		// C function pointer to receive log messages from the library
		pfnLoggerSink sink = nullptr;
		// Optional context parameter for the sink function; when consuming from C# you don't need that, pass IntPtr.Zero, delegates can capture things.
		void* context = nullptr;
		// Maximum log level to produce
		eLogLevel level;
	};
}  // namespace MeshRepair