#include "stdafx.h"
#include "Logger.h"
#include <cstdio>
using MeshRepair::eLogLevel;

inline bool Logger::shouldLog( eLogLevel lvl ) const
{
	return nullptr != setup.sink && (uint8_t)lvl <= (uint8_t)setup.level;
}

void Logger::logMessage( eLogLevel lvl, const char* pszFormat, std::va_list va ) const
{
	// Absolutely possible to support more than 128 characters per message, I just don't need to.
	constexpr size_t bufferSize = 128;
	char buf[ bufferSize ];
	int cc = std::vsnprintf( buf, bufferSize, pszFormat, va );
	if( cc < 0 )
		return;	 // fail
	if( cc >= bufferSize )
		buf[ bufferSize - 1 ] = '\0';
	else
		buf[ cc ] = '\0';
	setup.sink( setup.context, lvl, buf );
}

#define LOG_MESSAGE_IMPL( lvl )             \
	if( shouldLog( lvl ) )                  \
	{                                       \
		std::va_list args;                  \
		va_start( args, pszFormat );        \
		logMessage( lvl, pszFormat, args ); \
		va_end( args );                     \
	}

void Logger::logError( const char* pszFormat, ... ) const
{
	LOG_MESSAGE_IMPL( eLogLevel::Error );
}

void Logger::logWarning( const char* pszFormat, ... ) const
{
	LOG_MESSAGE_IMPL( eLogLevel::Warning );
}

void Logger::logInfo( const char* pszFormat, ... ) const
{
	LOG_MESSAGE_IMPL( eLogLevel::Info );
}

void Logger::logDebug( const char* pszFormat, ... ) const
{
	LOG_MESSAGE_IMPL( eLogLevel::Debug );
}

Logger::Logger( const MeshRepair::sLoggerSetup& sls )
	: setup( sls )
{
}