#include "stdafx.h"
#include <atlcoll.h>
#include <atlfile.h>
#include "IndexedMesh.h"
#include "simdUtils.h"

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

	struct HashKeyTraits : public CElementTraits<float3>
	{
		static inline bool CompareElements( const float3& a, const float3& b )
		{
			// We want 0.0 to compare unequal to -0.0, that's why integers
			__m128i v1 = loadInt3( &a );
			__m128i v2 = loadInt3( &b );
			__m128i xx = _mm_xor_si128( v1, v2 );
			return (bool)_mm_testz_si128( xx, xx );
		}

		static inline ULONG Hash( const float3& pos )
		{
			__m128i vec = loadInt3( &pos );
			// Set W lane to 1.0
			vec = _mm_insert_epi32( vec, 1, 3 );
			// https://stackoverflow.com/a/1646913/126995
			const __m128i mul = _mm_setr_epi32( 31, 31 * 31, 31 * 31 * 31, 31 * 31 * 31 * 31 );
			vec = _mm_mullo_epi32( vec, mul );
			vec = _mm_add_epi32( vec, _mm_srli_si128( vec, 8 ) );

			uint32_t y = (uint32_t)_mm_extract_epi32( vec, 1 );
			uint32_t x = (uint32_t)_mm_cvtsi128_si32( vec );
			return x + y;
		}
	};

	class HashMap : public CAtlMap<float3, uint32_t, HashKeyTraits>
	{
	  public:
		HashMap()
			: CAtlMap( 17u, 0.75f, 0.25f, 2.25f, 1024 )
		{
		}

		uint32_t add( __m128 pos )
		{
			std::array<int, 4> arr;
			_mm_storeu_si128( (__m128i*)arr.data(), _mm_castps_si128( pos ) );
			const float3& f3 = *(const float3*)( arr.data() );
			auto p = CAtlMap::Lookup( f3 );
			if( nullptr != p )
				return p->m_value;
			size_t idx = CAtlMap::GetCount();
			assert( idx < UINT_MAX );
			CAtlMap::SetAt( f3, (uint32_t)idx );
			return (uint32_t)idx;
		}

		HRESULT serialize( std::vector<float3>& vertices ) const
		{
			try
			{
				vertices.resize( CAtlMap::GetCount() );
			}
			catch( const std::bad_alloc& )
			{
				return E_OUTOFMEMORY;
			}

			float3* const rdi = vertices.data();
			for( POSITION pos = CAtlMap::GetStartPosition(); nullptr != pos; )
			{
				const auto* rsi = CAtlMap::GetNext( pos );
				rdi[ rsi->m_value ] = rsi->m_key;
			}
			return S_OK;
		}
	};
}  // namespace

HRESULT IndexedMesh::loadBinaryStl( LPCTSTR path )
{
	vertices.clear();
	triangles.clear();

	CAtlFile file;
	CHECK( file.Create( path, GENERIC_READ, FILE_SHARE_READ, OPEN_EXISTING ) );
	ULONGLONG fileSize = 0;
	CHECK( file.GetSize( fileSize ) );
	if( fileSize < 84 + sizeof( StlTriangle ) )
		return E_INVALIDARG;

	CAtlFileMapping<char> mapping;
	CHECK( mapping.MapFile( file ) );
	const char* ptr = mapping;
	if( 0 == StrCmpNIA( ptr, "solid", 5 ) )
		return E_INVALIDARG;  // We don't support text STL files

	ptr += 80;
	const uint32_t countTriangles = *(const uint32_t*)( ptr );
	if( fileSize != 84 + countTriangles * sizeof( StlTriangle ) )
		return E_INVALIDARG;

	const StlTriangle* rsi = (const StlTriangle*)( ptr + 4 );

	try
	{
		triangles.resize( countTriangles );
		std::array<uint32_t, 3>* rdi = triangles.data();
		HashMap hash;

		for( const StlTriangle* rsiEnd = rsi + countTriangles; rsi < rsiEnd; rsi++ )
		{
			const StlTriangle& tri = *rsi;
			const __m128 v0 = loadFloat3( &tri.vertices[ 0 ] );
			const __m128 v1 = loadFloat3( &tri.vertices[ 1 ] );
			const __m128 v2 = loadFloat3( &tri.vertices[ 2 ] );
			if( vectorEqual( v0, v1 ) || vectorEqual( v0, v2 ) || vectorEqual( v1, v2 ) )
				continue;
			std::array<uint32_t, 3>& dest = *rdi;
			dest[ 0 ] = hash.add( v0 );
			dest[ 1 ] = hash.add( v1 );
			dest[ 2 ] = hash.add( v2 );
			rdi++;
		}

		triangles.resize( rdi - triangles.data() );
		return hash.serialize( vertices );
	}
	catch( std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}
	catch( const CAtlException& ex )
	{
		return ex;
	}
}

namespace
{
	inline __m128 computeNormal( __m128 v0, __m128 v1, __m128 v2 )
	{
		__m128 e1 = _mm_sub_ps( v1, v0 );
		__m128 e2 = _mm_sub_ps( v2, v0 );
		__m128 cp = DirectX::XMVector3Cross( e1, e2 );
		return DirectX::XMVector3NormalizeEst( cp );
	}

	inline void makeTriangle( StlTriangle& rdi, const float3* vb, const std::array<uint32_t, 3>& tri )
	{
		__m128 v0 = loadFloat3( vb + tri[ 0 ] );
		__m128 v1 = loadFloat3( vb + tri[ 1 ] );
		__m128 v2 = loadFloat3( vb + tri[ 2 ] );
		__m128 norm = computeNormal( v0, v1, v2 );
		storeFloat3_x4( norm, v0, v1, v2, &rdi.normal );
	}
}  // namespace

HRESULT IndexedMesh::saveBinaryStl( LPCTSTR path ) const
{
	if( triangles.empty() )
		return OLE_E_BLANK;
	if( triangles.size() > UINT_MAX )
		return DISP_E_OVERFLOW;

	CAtlFile file;
	CHECK( file.Create( path, GENERIC_READ | GENERIC_WRITE, 0, CREATE_ALWAYS ) );
	CHECK( file.SetSize( 84 + triangles.size() * sizeof( StlTriangle ) ) );

	CAtlFileMapping<char> mapping;
	CHECK( mapping.MapFile( file, 0, 0, PAGE_READWRITE, FILE_MAP_ALL_ACCESS ) );

	char* ptr = (char*)mapping.GetData();
	*(uint32_t*)( ptr + 80 ) = (uint32_t)triangles.size();

	const float3* const vb = vertices.data();
	StlTriangle* rdi = (StlTriangle*)( ptr + 84 );

	for( const auto& tri : triangles )
	{
		makeTriangle( *rdi, vb, tri );
		rdi++;
	}

	CHECK( mapping.Unmap() );
	CHECK( file.Flush() );
	file.Close();
	return S_OK;
}