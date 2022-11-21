#include "stdafx.h"
#include "ResultMesh.h"
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
	const size_t countTriangles = F_sf.rows();
	memcpy( rdi, F_sf.data(), countTriangles * 12 );
	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh::getVertexBufferFP32( float* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );

	const size_t countVertices = V_sf.rows();
	const size_t countValues = countVertices * 3;

	const double* rsi = V_sf.data();
	const double* const rsiEnd = rsi + countValues;
	const double* const rsiEndAligned = rsi + ( countValues / 4 ) * 4;

	for( ; rsi < rsiEndAligned; rsi += 4, rdi += 4 )
	{
		__m256d vd = _mm256_loadu_pd( rsi );
		__m128 vf = _mm256_cvtpd_ps( vd );
		_mm_storeu_ps( rdi, vf );
	}

#pragma loop( no_vector )
	for( ; rsi < rsiEnd; rsi++, rdi++ )
	{
		double d = *rsi;
		float f = (float)d;
		*rdi = f;
	}

	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh::getVertexBufferFP64( double* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );

	const size_t countVertices = V_sf.rows();
	const size_t countValues = countVertices * 3;
	memcpy( rdi, V_sf.data(), countValues * 8 );

	return S_OK;
}