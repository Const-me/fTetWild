#include "stdafx.h"
#include "TriMeshWrapper.h"
#include "../miscUtils.h"

HRESULT TriMeshWrapper::assignVertices( size_t count, const float* vb )
{
	const size_t totalFloats = count * 3;
	GEO::vector<double> vec;
	try
	{
		vec.resize( totalFloats );
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}

	upcastFloats( vec.data(), totalFloats, vb );
	m.vertices.assign_points( vec, 3, true );
	return S_OK;
}

HRESULT TriMeshWrapper::assignTriangles( size_t count, const uint32_t* ib )
{
	const size_t totalIntegers = count * 3;
	GEO::vector<GEO::index_t> vec;
	try
	{
		vec.resize( totalIntegers );
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}

	static_assert( sizeof( GEO::index_t ) == 4, "GEO2::index_t expected to be 4 bytes" );
	memcpy( vec.data(), ib, totalIntegers * 4 );
	m.facets.assign_triangle_mesh( vec, true );
	return S_OK;
}

HRESULT TriMeshWrapper::copyData( std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib ) const
{
	// Extract the data, store in different types
	try
	{
		vb.resize( m.vertices.nb() );
		ib.resize( m.facets.nb() );
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}

	for( size_t i = 0; i < vb.size(); i++ )
	{
		const GEO::vec3& src = m.vertices.point( (GEO::index_t)i );
		copyDouble3( vb[ i ].data(), &src.x );
	}

	for( size_t i = 0; i < ib.size(); i++ )
	{
		uint32_t* rdi = (uint32_t*)&ib[ i ];
		rdi[ 0 ] = m.facets.vertex( (GEO::index_t)i, 0 );
		rdi[ 1 ] = m.facets.vertex( (GEO::index_t)i, 1 );
		rdi[ 2 ] = m.facets.vertex( (GEO::index_t)i, 2 );
	}

	return S_OK;
}