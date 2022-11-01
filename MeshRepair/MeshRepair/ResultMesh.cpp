#include "stdafx.h"
#include "ResultMesh.h"
using namespace MeshRepair;

HRESULT ResultMesh::ensureGoodMesh() const
{
	if( V_sf.cols() != 3 || F_sf.cols() != 3 || V_sf.rows() == 0 || F_sf.rows() == 0 )
		return OLE_E_BLANK;

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

	const int* px = &F_sf( 0, 0 );
	const int* py = &F_sf( 0, 1 );
	const int* pz = &F_sf( 0, 2 );
	const int* const pxEnd = px + ( F_sf.rows() - 1 );

	for( ; px < pxEnd; px++, py++, pz++, rdi += 3 )
	{
		// Most modern CPUs can do twice as many loads per cycle, compared to stores.
		// Optimizing for smaller count of stores.
		__m128i v = _mm_loadu_si32( px );
		v = _mm_insert_epi32( v, *py, 1 );
		v = _mm_insert_epi32( v, *pz, 2 );
		_mm_storeu_si128( (__m128i*)rdi, v );
	}

	// For the last element we can't store 16 unaligned bytes instead of 12, would corrupt memory.
	{
		rdi[ 0 ] = *(const uint32_t*)px;
		rdi[ 1 ] = *(const uint32_t*)py;
		rdi[ 2 ] = *(const uint32_t*)pz;
	}
	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh::getVertexBufferFP32( float* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );

	const double* px = &V_sf( 0, 0 );
	const double* py = &V_sf( 0, 1 );
	const double* pz = &V_sf( 0, 2 );
	const double* const pxEnd = px + ( V_sf.rows() - 1 );

	for( ; px < pxEnd; px++, py++, pz++, rdi += 3 )
	{
		__m128d d = _mm_load_sd( px );
		d = _mm_loadh_pd( d, py );
		__m128 f = _mm_cvtpd_ps( d );

		d = _mm_load_sd( pz );
		__m128 fz = _mm_cvtpd_ps( d );
		f = _mm_movelh_ps( f, fz );
		_mm_storeu_ps( rdi, f );
	}

	// For the last element we can't store 16 unaligned bytes instead of 12, would corrupt memory.
	{
		__m128d d = _mm_load_sd( px );
		d = _mm_loadh_pd( d, py );
		__m128 f = _mm_cvtpd_ps( d );
		_mm_store_sd( (double*)rdi, _mm_castps_pd( f ) );

		d = _mm_load_sd( pz );
		f = _mm_cvtpd_ps( d );
		_mm_store_ss( rdi + 2, f );
	}
	return S_OK;
}

HRESULT COMLIGHTCALL ResultMesh::getVertexBufferFP64( double* rdi ) const noexcept
{
	CHECK( ensureGoodMesh() );

	const double* px = &V_sf( 0, 0 );
	const double* py = &V_sf( 0, 1 );
	const double* pz = &V_sf( 0, 2 );
	const double* const pxEnd = px + V_sf.rows();

	for( ; px < pxEnd; px++, py++, pz++, rdi += 3 )
	{
		__m128d d = _mm_load_sd( px );
		d = _mm_loadh_pd( d, py );
		_mm_storeu_pd( rdi, d );
		rdi[ 2 ] = *pz;
	}

	return S_OK;
}