#include "stdafx.h"
#include "ConsoleLogSink.h"

namespace
{
	using MeshRepair::eLogLevel;
	// https://github.com/Const-me/vis_avs_dx/blob/master/avs_dx/DxVisuals/Interop/ConsoleLogger.cpp

	constexpr uint16_t defaultAttributes = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE;

	inline uint16_t textAttributes( eLogLevel lvl )
	{
		switch( lvl )
		{
		case eLogLevel::Error:
			return FOREGROUND_RED | FOREGROUND_INTENSITY;
		case eLogLevel::Warning:
			return FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY;
		case eLogLevel::Info:
			return FOREGROUND_GREEN | FOREGROUND_INTENSITY;
		case eLogLevel::Debug:
			return FOREGROUND_BLUE | FOREGROUND_INTENSITY;
		}
		return defaultAttributes;
	}

	static void setTextAttributes( uint16_t ta )
	{
		const HANDLE h = GetStdHandle( STD_OUTPUT_HANDLE );
		if( nullptr != h )
			SetConsoleTextAttribute( h, ta );
	}
}  // namespace

ConsoleLogSink::SetupColors::SetupColors( eLogLevel level )
{
	setTextAttributes( textAttributes( level ) );
}

ConsoleLogSink::SetupColors::~SetupColors()
{
	setTextAttributes( defaultAttributes );
}