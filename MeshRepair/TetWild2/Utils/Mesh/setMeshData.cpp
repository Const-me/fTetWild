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
	return mesh.assignVertices( count, vb );
}

HRESULT assignMeshTriangles( GEO2::Mesh& mesh, size_t count, const uint32_t* ib ) 
{
	return mesh.assignTriangles( count, ib );
}

HRESULT copyMeshData( const GEO2::Mesh& mesh, std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib ) 
{
	return mesh.copyData( vb, ib );
}