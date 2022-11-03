#pragma once
#include <geogram/mesh/mesh_AABB.h>

class TetraMeshWrapper : private GEO::Mesh
{
  public:
	GEO::Mesh& unwrap()
	{
		return *this;
	}

	template<class Lambda>
	void generateVertices( uint32_t count, Lambda lambda )
	{
		vertices.clear();
		vertices.create_vertices( count );
		for( uint32_t i = 0; i < count; i++ )
			vertices.point( i ) = lambda( i );
	}

	template<class Lambda>
	void generateElements( uint32_t count, Lambda lambda )
	{
		cells.create_tets( (int)count );
		for( uint32_t i = 0; i < count; i++ )
		{
			const __m128i vec = lambda( i );
			cells.set_vertex( i, 0, _mm_cvtsi128_si32( vec ) );
			cells.set_vertex( i, 1, _mm_extract_epi32( vec, 1 ) );
			cells.set_vertex( i, 2, _mm_extract_epi32( vec, 2 ) );
			cells.set_vertex( i, 3, _mm_extract_epi32( vec, 3 ) );
		}
	}
};

class MeshCellsAABBWrapper : private GEO::MeshCellsAABB
{
  public:
	MeshCellsAABBWrapper( TetraMeshWrapper& mesh, bool reorder = true )
		: GEO::MeshCellsAABB( mesh.unwrap(), reorder )
	{
	}

	uint32_t containingElement( const GEO::vec3& p, bool exact = true ) const
	{
		return MeshCellsAABB::containing_tet( p, exact );
	}

	static constexpr uint32_t NO_TET = UINT_MAX;
};