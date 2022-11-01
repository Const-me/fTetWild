#include "stdafx.h"
#include "SourceMesh.h"
using namespace MeshRepair;

namespace
{
	// Set vertex buffer of the mesh, upcasting coordinates to FP64
	HRESULT assignMeshVertices( GEO::Mesh& mesh, size_t count, const float* vb )
	{
		const size_t totalFloats = count * 3;
		GEO::vector<double> vec;
		try
		{
			vec.resize( totalFloats );
		}
		catch( const std::bad_alloc& ba )
		{
			return E_OUTOFMEMORY;
		}

		const float* rsiEnd = vb + totalFloats;
		double* rdi = vec.data();
		const float* rsiEndAligned = vb + ( totalFloats / 4 ) * 4;
		for( ; vb < rsiEndAligned; rdi += 4 )
		{
			__m128 v = _mm_loadu_ps( vb );
			vb += 4;
#ifdef __AVX__
			__m256d d = _mm256_cvtps_pd( v );
			_mm256_storeu_pd( rdi, d );
#else
			__m128d d = _mm_cvtps_pd( v );
			v = _mm_unpackhi_ps( v, v );
			_mm_storeu_pd( rdi, d );
			d = _mm_cvtps_pd( v );
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

		mesh.vertices.assign_points( vec, 3, true );
		return S_OK;
	}

	// Set index buffer of the mesh
	HRESULT assignMeshTriangles( GEO::Mesh& mesh, size_t count, const uint32_t* ib ) 
	{
		const size_t totalIntegers = count * 3;
		GEO::vector<GEO::index_t> vec;
		try
		{
			vec.resize( totalIntegers );
		}
		catch( const std::bad_alloc& ba )
		{
			return E_OUTOFMEMORY;
		}

		static_assert( sizeof( GEO::index_t ) == 4, "GEO::index_t expected to be 4 bytes" );
		memcpy( vec.data(), ib, count * 4 );
		mesh.facets.assign_triangle_mesh( vec, true );
		return S_OK;
	}
}  // namespace

HRESULT SourceMesh::createMesh( uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib )
{
	if( countVertices == 0 || countTriangles == 0 )
		return E_INVALIDARG;
	if( nullptr == vb || nullptr == ib )
		return E_POINTER;

	CHECK( assignMeshVertices( mesh, countVertices, vb ) );
	CHECK( assignMeshTriangles( mesh, countTriangles, ib ) );

	return S_OK;
}