#include "stdafx.h"
#include "TriangleMesh.h"
#include "setMeshData.h"

// Set vertex buffer of the mesh, upcasting coordinates to FP64
HRESULT TriangleMesh::assignVertices( size_t count, const float* vb )
{
	try
	{
		vertices.resize( count );
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}
	upcastFloats( (double*)vertices.data(), count * 3, vb );
	return S_OK;
}

HRESULT TriangleMesh::assignTriangles( size_t count, const uint32_t* ib )
{
	static_assert( sizeof( GEO::vec3i ) == 12, "sizeof" );

	try
	{
		const GEO::vec3i* rsi = (const GEO::vec3i*)ib;
		triangles.assign( rsi, rsi + count );
		return S_OK;
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}
}

double TriangleMesh::boxDiagonal() const
{
	__m128d xyMin = _mm_set1_pd( DBL_MAX );
	__m128d zMin = xyMin;
	__m128d xyMax = _mm_sub_pd( _mm_setzero_pd(), xyMin );
	__m128d zMax = xyMax;

	for( const GEO::vec3& vert : vertices )
	{
		const double* rsi = (const double*)&vert;
		__m128d v = _mm_loadu_pd( rsi );
		xyMin = _mm_min_pd( xyMin, v );
		xyMax = _mm_max_pd( xyMax, v );

		v = _mm_load_sd( rsi + 2 );
		zMin = _mm_min_pd( zMin, v );
		zMax = _mm_max_pd( zMax, v );
	}

	__m128d xy = _mm_sub_pd( xyMax, xyMin );
	__m128d z = _mm_sub_pd( zMax, zMin );

	xy = _mm_dp_pd( xy, xy, 0b00110001 );
	z = _mm_mul_sd( z, z );

	__m128d res = _mm_add_sd( xy, z );
	res = _mm_sqrt_sd( res, res );
	return _mm_cvtsd_f64( res );
}

void TriangleMesh::clearMesh( bool keepMemory )
{
	vertices.clear();
	triangles.clear();
	if( keepMemory )
		return;

	vertices.shrink_to_fit();
	triangles.shrink_to_fit();
}