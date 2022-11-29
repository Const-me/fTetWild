#pragma once
#include <cstdarg>
#include "../../MeshRepair/API/loggerApi.h"

class Logger
{
	MeshRepair::sLoggerSetup setup;

	inline bool shouldLog( MeshRepair::eLogLevel lvl ) const;
	void logMessage( MeshRepair::eLogLevel lvl, const char* pszFormat, std::va_list va ) const;

  public:
	Logger( const MeshRepair::sLoggerSetup& sls );

	void logError( const char* pszFormat, ... ) const;
	void logWarning( const char* pszFormat, ... ) const;
	void logInfo( const char* pszFormat, ... ) const;
	void logDebug( const char* pszFormat, ... ) const;
};