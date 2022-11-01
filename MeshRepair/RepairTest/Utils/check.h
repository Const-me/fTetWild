#pragma once

#define CHECK( hr )                  \
	{                                \
		const HRESULT __hr = ( hr ); \
		if ( FAILED( __hr ) )        \
			return __hr;             \
	}