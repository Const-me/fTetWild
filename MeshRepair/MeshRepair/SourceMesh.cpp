#include "stdafx.h"
#include "SourceMesh.h"
#include <geogram/mesh/mesh_reorder.h>
using namespace MeshRepair;

namespace
{
	static void upcastFloats( double* rdi, size_t length, const float* vb ) 
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

		upcastFloats( vec.data(), totalFloats, vb );
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
		memcpy( vec.data(), ib, totalIntegers * 4 );
		mesh.facets.assign_triangle_mesh( vec, true );
		return S_OK;
	}

	inline void copyDouble3( double* rdi, const double* rsi )
	{
		__m128d v = _mm_loadu_pd( rsi );
		_mm_storeu_pd( rdi, v );
		rdi[ 2 ] = rsi[ 2 ];
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

	// See MeshIO::load_mesh function in the original code
	GEO::mesh_reorder( mesh, GEO::MESH_ORDER_MORTON );

	// Extract the data, store in different types
	try
	{
		input_vertices.resize( mesh.vertices.nb() );
		input_faces.resize( mesh.facets.nb() );
		input_tags.clear();
		input_tags.resize( input_faces.size(), 0 );
	}
	catch( const std::bad_alloc& ba )
	{
		return E_OUTOFMEMORY;
	}

    for( size_t i = 0; i < input_vertices.size(); i++ ) 
	{
		const GEO::vec3& src = mesh.vertices.point( i );
		copyDouble3( input_vertices[ i ].data(), &src.x );
	}

	for( size_t i = 0; i < input_faces.size(); i++ ) 
	{
		uint32_t* rdi = (uint32_t*)&input_faces[ i ];
		rdi[ 0 ] = mesh.facets.vertex( i, 0 );
		rdi[ 1 ] = mesh.facets.vertex( i, 1 );
		rdi[ 2 ] = mesh.facets.vertex( i, 2 );
	}

	return S_OK;
}