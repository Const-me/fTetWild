#include "stdafx.h"
#include "writeStl.h"
#include <DirectXMath.h>
#include <atlfile.h>

namespace
{
	using float3 = DirectX::XMFLOAT3;

#pragma pack( push )
#pragma pack( 1 )
	struct StlTriangle
	{
		float3 normal;
		std::array<float3, 3> vertices;
		uint16_t unused;
	};
#pragma pack( pop )

	inline void throwIfFailed( HRESULT hr )
	{
		if( SUCCEEDED( hr ) )
			return;
		__debugbreak();
	}

	inline __m128 loadFloat3( const floatTetWild::Vector3* rsi )
	{
		const double* p = (const double*)rsi;
		__m128d d = _mm_loadu_pd( p );
		__m128 f = _mm_cvtpd_ps( d );
		d = _mm_load_sd( p + 2 );
		return _mm_movelh_ps( f, _mm_cvtpd_ps( d ) );
	}

	inline __m128 computeNormal( __m128 v0, __m128 v1, __m128 v2 )
	{
		__m128 e1 = _mm_sub_ps( v1, v0 );
		__m128 e2 = _mm_sub_ps( v2, v0 );
		__m128 cp = DirectX::XMVector3Cross( e1, e2 );
		return DirectX::XMVector3NormalizeEst( cp );
	}

	// Make a mask for _mm_insert_ps instruction
	constexpr int insertMask( int source, int dest, int zero = 0 )
	{
		assert( source >= 0 && source < 4 );
		assert( dest >= 0 && dest < 4 );
		assert( zero >= 0 && zero < 16 );
		return ( source << 6 ) | ( dest << 4 ) | zero;
	}

	// Store 12 floats from XYZ lanes of 4 vectors, with 3 store instructions.
	// The function writes 12 floats = 48 bytes to the pointer.
	inline void storeFloat3_x4( __m128 a, __m128 b, __m128 c, __m128 d, void* dest )
	{
		float* const rdi = (float*)dest;

		// a.xyz, b.x
		__m128 tmp = _mm_insert_ps( a, b, insertMask( 0, 3, 0 ) );
		_mm_storeu_ps( rdi, tmp );

		// b.yz, c.xy
		tmp = _mm_shuffle_ps( b, c, _MM_SHUFFLE( 1, 0, 2, 1 ) );
		_mm_storeu_ps( rdi + 4, tmp );

		// c.z, d.xyz
		tmp = _mm_permute_ps( d, _MM_SHUFFLE( 2, 1, 0, 0 ) );  // Permute into d.xxyz
		tmp = _mm_insert_ps( tmp, c, insertMask( 2, 0, 0 ) );
		_mm_storeu_ps( rdi + 8, tmp );
	}

	inline void makeTriangle( StlTriangle& rdi, const floatTetWild::Vector3* vb, const floatTetWild::Vector3i& tri )
	{
		__m128 v0 = loadFloat3( vb + tri[ 0 ] );
		__m128 v1 = loadFloat3( vb + tri[ 1 ] );
		__m128 v2 = loadFloat3( vb + tri[ 2 ] );
		__m128 norm = computeNormal( v0, v1, v2 );
		storeFloat3_x4( norm, v0, v1, v2, &rdi.normal );
	}
}  // namespace

HRESULT writeStl( const std::vector<floatTetWild::Vector3> vertices, const std::vector<floatTetWild::Vector3i>& triangles, LPCTSTR path )
{
	if( triangles.empty() )
		return OLE_E_BLANK;
	if( triangles.size() > UINT_MAX )
		return DISP_E_OVERFLOW;

	CAtlFile file;
	throwIfFailed( file.Create( path, GENERIC_READ | GENERIC_WRITE, 0, CREATE_ALWAYS ) );
	throwIfFailed( file.SetSize( 84 + triangles.size() * sizeof( StlTriangle ) ) );

	CAtlFileMapping<char> mapping;
	throwIfFailed( mapping.MapFile( file, 0, 0, PAGE_READWRITE, FILE_MAP_ALL_ACCESS ) );

	char* ptr = (char*)mapping.GetData();
	*(uint32_t*)( ptr + 80 ) = (uint32_t)triangles.size();

	const floatTetWild::Vector3* const vb = vertices.data();
	StlTriangle* rdi = (StlTriangle*)( ptr + 84 );

	for( const auto& tri : triangles )
	{
		makeTriangle( *rdi, vb, tri );
		rdi++;
	}

	throwIfFailed( mapping.Unmap() );
	throwIfFailed( file.Flush() );
	file.Close();
	return S_OK;
}