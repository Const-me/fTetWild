#include "stdafx.h"
#include "ConsoleLogSink.h"
#include <stdio.h>

inline void ConsoleLogSink::logMessage( MeshRepair::eLogLevel lvl, const char* message )
{
	std::lock_guard<std::mutex> lk { mutex };
	SetupColors col { lvl };
	printf( "%s\n", message );
}

void __cdecl ConsoleLogSink::logMessageStatic( void* context, MeshRepair::eLogLevel lvl, const char* message )
{
	ConsoleLogSink* cls = (ConsoleLogSink*)context;
	cls->logMessage( lvl, message );
}

ConsoleLogSink::ConsoleLogSink( MeshRepair::eLogLevel level )
{
	setup.sink = &logMessageStatic;
	setup.context = this;
	setup.level = level;
}