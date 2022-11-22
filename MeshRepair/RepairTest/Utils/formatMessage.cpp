#include "stdafx.h"
#include "formatMessage.h"

namespace
{
	constexpr DWORD flags = FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS;
	constexpr DWORD langId = MAKELANGID( LANG_NEUTRAL, SUBLANG_DEFAULT );

	inline DWORD formatMessageW( HRESULT hr, wchar_t** ppStr )
	{
		return ::FormatMessageW( flags, nullptr, hr, langId, (LPWSTR)ppStr, 0, nullptr );
	}
}

CString formatMessage( HRESULT hr )
{
	wchar_t* lpMsgBuf = nullptr;
	const DWORD fmLength = formatMessageW( hr, &lpMsgBuf );
	if( 0 != fmLength )
	{
		CString msg { lpMsgBuf, (int)fmLength };
		LocalFree( lpMsgBuf );
		msg.Trim( L"\r\n\t " );
		return msg;
	}
	LocalFree( lpMsgBuf );

	// Fallback message
	CStringW unknown;
	unknown.Format( L"Unknown error code %i (0x%08X)", (int)hr, (int)hr );
	return unknown;
}