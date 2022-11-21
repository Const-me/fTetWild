#include "stdafx.h"
#include "ResultMesh.h"
#include "Utils/downcastMesh.h"
using namespace MeshRepair;

HRESULT ResultMesh::ensureGoodMesh() const
{
	if( V_sf.rows() == 0 || F_sf.rows() == 0 )
		return OLE_E_BLANK;
	if( 3 != V_sf.rowStride() || 3 != F_sf.rowStride() )
		return E_UNEXPECTED;

	// If you have meshes with >2G vertices or triangles, refactor the library API replacing sizes with uint64_t
	if( V_sf.rows() > INT_MAX || F_sf.rows() > INT_MAX )
		return DISP_E_OVERFLOW;

	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh::getSize( uint32_t& vertices, uint32_t& triangles ) const noexcept
{
	CHECK( ensureGoodMesh() );
	vertices = (uint32_t)V_sf.rows();
	triangles = (uint32_t)F_sf.rows();
	return S_OK;
}

// By default, Eigen matrices are column-major, transposing numbers on the fly.

HRESULT COMLIGHTCALL ResultMesh::getIndexBuffer( uint32_t* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );
	if( nullptr == rdi )
		return E_POINTER;
	const size_t countTriangles = F_sf.rows();
	memcpy( rdi, F_sf.data(), countTriangles * 12 );
	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh::getVertexBufferFP32( float* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );
	if( nullptr == rdi )
		return E_POINTER;

	const size_t countVertices = V_sf.rows();
	const size_t countValues = countVertices * 3;
	downcastFloats( rdi, countValues, V_sf.data() );
	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh::getVertexBufferFP64( double* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );
	if( nullptr == rdi )
		return E_POINTER;

	const size_t countVertices = V_sf.rows();
	const size_t countValues = countVertices * 3;
	memcpy( rdi, V_sf.data(), countValues * 8 );

	return S_OK;
}

HRESULT ResultMesh32::ensureGoodMesh() const
{
	if( vertexBuffer.empty() || indexBuffer.empty() )
		return OLE_E_BLANK;
	if( vertexBuffer.size() > INT_MAX || indexBuffer.size() > INT_MAX )
		return DISP_E_OVERFLOW;
	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh32::getSize( uint32_t& vertices, uint32_t& triangles ) const noexcept
{
	CHECK( ensureGoodMesh() );
	vertices = (uint32_t)vertexBuffer.size();
	triangles = (uint32_t)indexBuffer.size();
	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh32::getIndexBuffer( uint32_t* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );
	if( nullptr == rdi )
		return E_POINTER;
	const size_t countTriangles = indexBuffer.size();
	memcpy( rdi, indexBuffer.data(), countTriangles * 12 );
	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh32::getVertexBufferFP32( float* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );
	if( nullptr == rdi )
		return E_POINTER;

	const size_t countVertices = vertexBuffer.size();
	memcpy( rdi, vertexBuffer.data(), countVertices * 12 );
	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh32::getVertexBufferFP64( double* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );
	if( nullptr == rdi )
		return E_POINTER;

	const size_t countVertices = vertexBuffer.size();
	const size_t countValues = countVertices * 3;

	const float* rsi = (float*)vertexBuffer.data();
	const float* const rsiEnd = rsi + countValues;
	const float* const rsiEndAligned = rsi + ( countValues / 4 ) * 4;

	for( ; rsi < rsiEndAligned; rsi += 4, rdi += 4 )
	{
		__m128 vf = _mm_loadu_ps( rsi );
		__m256d vd = _mm256_cvtps_pd( vf );
		_mm256_storeu_pd( rdi, vd );
	}

#pragma loop( no_vector )
	for( ; rsi < rsiEnd; rsi++, rdi++ )
		*rdi = *rsi;

	return S_OK;
}