#pragma once
#include "MeshBase.h"

// Triangular surface mesh in 3D
class TriangleMesh : public MeshBase
{
	std::vector<GEO2::vec3i> triangles;

  public:
	using vec3i = GEO2::vec3i;

	// Compute 3.0 * specified coordinate of the center of the triangle
	template<int COORD>
	inline double triangleCenterX3( uint32_t idxTri ) const;

	// Set index buffer of the mesh
	HRESULT assignTriangles( size_t count, const uint32_t* ib );
	// Copy vertex and index buffer out of the mesh
	HRESULT copyData( std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib ) const;

	uint32_t countTriangles() const
	{
		return triangles.size();
	}

	const vec3i* trianglePointer() const
	{
		return triangles.data();
	}
	vec3i* trianglePointer()
	{
		return triangles.data();
	}
	const vec3i& getTriangle( uint32_t t ) const
	{
		return triangles[ t ];
	}

	void reorderMorton();

	// Compute diagonal of the bounding box
	double boxDiagonal() const;

	void clearMesh( bool keepMemory = true );

	void createTriangles( size_t count )
	{
		triangles.resize( count );
	}

	// Get constant pointers to all vertices of the specified triangle
	void getTriangleVertices( uint32_t tri, const vec3** p1, const vec3** p2, const vec3** p3 ) const
	{
		const vec3i& t = getTriangle( tri );
		*p1 = &getVertex( t.x );
		*p2 = &getVertex( t.y );
		*p3 = &getVertex( t.z );
	}

	const vec3& getFirstTriangleVertex( uint32_t tri ) const
	{
		const vec3i& t = getTriangle( tri );
		return getVertex( t.x );
	}

	template<class Lambda>
	void generateTriangles( uint32_t count, Lambda lambda )
	{
		triangles.resize( count );
		vec3i* rdi = triangles.data();
		for( uint32_t i = 0; i < count; i++, rdi++ )
			*rdi = lambda( i );
	}
};