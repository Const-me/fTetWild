#pragma once
#include <atlstr.h>

CString formatMessage( HRESULT hr );

inline const wchar_t* cstr( const CString& s )
{
	return s;
}