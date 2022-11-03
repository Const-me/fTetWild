#include "stdafx.h"
#include "designTimeUtils.h"
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <atlstr.h>

namespace
{
	// Print float as hexadecimal literal compatible with C++ compiler
	static void printDouble( CStringA& rdi, double f )
	{
		if( f == 0.0 )
		{
			rdi.Append( std::signbit( f ) ? "-0.0" : "0.0" );
			return;
		}

		const int ol = rdi.GetLength();
		rdi.AppendFormat( "%a", f );

		// For optimal readability, trim trailing zeros from mantissa. This truncates "0x1.9999a80000000p-4" into "0x1.9999a8p-4"
		int idx = rdi.Find( "0p", ol );
		if( idx >= 0 )
		{
			const char* const ps = rdi;
			assert( ps[ idx ] == '0' );
			// end is index of the first character to keep
			const int end = idx + 1;

			while( true )
			{
				if( idx <= ol )
					break;
				int ip = idx - 1;
				if( ps[ ip ] != '0' )
					break;
				idx = ip;
			}
			// idx is index of the first character to overwrite
			assert( ps[ idx ] == '0' );
			assert( ps[ idx - 1 ] != '0' );

			const int nl = rdi.GetLength();
			char* p = rdi.GetBuffer();
			memmove( p + idx, p + end, nl - end );
			const int trimmedCharacters = end - idx;
			rdi.ReleaseBufferSetLength( nl - trimmedCharacters );
		}
	}
}

void dbgPrintCppLiteral( const char* name, double val )
{
	CStringA message;
	message.AppendFormat( "static const double %s = ", name );
	printDouble( message, val );
	message.Append( ";\r\n" );
	OutputDebugStringA( message );
}