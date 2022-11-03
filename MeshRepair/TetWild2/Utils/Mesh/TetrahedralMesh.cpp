#include "stdafx.h"
#include "TetrahedralMesh.h"

void TetrahedralMesh::assignVertices( size_t count, const double* vb )
{
	static_assert( alignof( GEO::vec3 ) == alignof( double ), "misalignmenet" );
	const GEO::vec3* rsi = (const GEO::vec3*)vb;
	vertices.assign( rsi, rsi + count );
}

void TetrahedralMesh::assignElements( size_t count, const uint32_t* ib )
{
	// vector.assign, may cause runtime crash in SSE builds, because source data is not guaranteed to be aligned by 16 bytes.
	// Thet's why memcpy() instead
	elements.resize( count );
	memcpy( elements.data(), ib, count * 16 );
}