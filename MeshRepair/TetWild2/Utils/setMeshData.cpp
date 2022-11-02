#include "stdafx.h"
#include "setMeshData.h"

void upcastFloats( double* rdi, size_t length, const float* vb )
{
	const float* rsiEnd = vb + length;
	const float* rsiEndAligned = vb + ( length / 4 ) * 4;
	for( ; vb < rsiEndAligned; rdi += 4 )
	{
		__m128 v = _mm_loadu_ps( vb );
		vb += 4;
#ifdef __AVX__
		__m256d d = _mm256_cvtps_pd( v );
		_mm256_storeu_pd( rdi, d );
#else
		// Upcast v.xy to FP64
		__m128d d = _mm_cvtps_pd( v );
		// Permute the source FP32 vector into v.zwzw
		v = _mm_movehl_ps( v, v );
		// Store first two FP64 numbers
		_mm_storeu_pd( rdi, d );
		// Upcast v.xy to FP64
		d = _mm_cvtps_pd( v );
		// Store remaining two FP64 numbers
		_mm_storeu_pd( rdi + 2, d );
#endif
	}

#pragma loop( no_vector )
	for( ; vb < rsiEnd; rdi++ )
	{
		float f = *vb;
		vb++;
		*rdi = f;
	}
}

HRESULT assignMeshVertices( GEO2::Mesh& mesh, size_t count, const float* vb ) 
{
	const size_t totalFloats = count * 3;
	GEO2::vector<double> vec;
	try
	{
		vec.resize( totalFloats );
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}

	upcastFloats( vec.data(), totalFloats, vb );
	mesh.vertices.assign_points( vec, 3, true );
	return S_OK;
}

HRESULT assignMeshTriangles( GEO2::Mesh& mesh, size_t count, const uint32_t* ib ) 
{
	const size_t totalIntegers = count * 3;
	GEO2::vector<GEO2::index_t> vec;
	try
	{
		vec.resize( totalIntegers );
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}

	static_assert( sizeof( GEO2::index_t ) == 4, "GEO2::index_t expected to be 4 bytes" );
	memcpy( vec.data(), ib, totalIntegers * 4 );
	mesh.facets.assign_triangle_mesh( vec, true );

	return S_OK;
}

namespace
{
	inline void copyDouble3( double* rdi, const double* rsi )
	{
		__m128d v = _mm_loadu_pd( rsi );
		_mm_storeu_pd( rdi, v );
		rdi[ 2 ] = rsi[ 2 ];
	}
}

HRESULT copyMeshData( const GEO2::Mesh& mesh, std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib ) 
{
	// Extract the data, store in different types
	try
	{
		vb.resize( mesh.vertices.nb() );
		ib.resize( mesh.facets.nb() );
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}

	for( size_t i = 0; i < vb.size(); i++ )
	{
		const GEO2::vec3& src = mesh.vertices.point( (GEO2::index_t)i );
		copyDouble3( vb[ i ].data(), &src.x );
	}

	for( size_t i = 0; i < ib.size(); i++ )
	{
		uint32_t* rdi = (uint32_t*)&ib[ i ];
		rdi[ 0 ] = mesh.facets.vertex( (GEO2::index_t)i, 0 );
		rdi[ 1 ] = mesh.facets.vertex( (GEO2::index_t)i, 1 );
		rdi[ 2 ] = mesh.facets.vertex( (GEO2::index_t)i, 2 );
	}

	return S_OK;
}