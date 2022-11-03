#pragma once
#include <geogram/basic/geometry.h>
#include "../../src/Types.hpp"

class TriangleMesh
{
	GEO::vector<GEO::vec3> vertices;
	GEO::vector<GEO::vec3i> triangles;

  public:
	// Compute 3.0 * specified coordinate of the center of the triangle
	template<int COORD>
	inline double triangleCenterX3( uint32_t idxTri ) const;

	// Set vertex buffer of the mesh, upcasting coordinates to FP64
	HRESULT assignVertices( size_t count, const float* vb );
	// Set index buffer of the mesh
	HRESULT assignTriangles( size_t count, const uint32_t* ib );
	// Copy vertex and index buffer out of the mesh
	HRESULT copyData( std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib ) const;

	uint32_t countVertices() const
	{
		return vertices.size();
	}

	uint32_t countTriangles() const
	{
		return triangles.size();
	}

	const GEO::vec3* vertexPointer() const
	{
		return vertices.data();
	}
	const GEO::vec3i* trianglePointer() const
	{
		return triangles.data();
	}
	GEO::vec3* vertexPointer()
	{
		return vertices.data();
	}
	GEO::vec3i* trianglePointer()
	{
		return triangles.data();
	}
	const GEO::vec3i& getTriangle( uint32_t t ) const
	{
		return triangles[ t ];
	}
	const GEO::vec3& getVertex( uint32_t v ) const
	{
		return vertices[ v ];
	}

	void reorderMorton();

	// Compute diagonal of the bounding box
	double boxDiagonal() const;

	void clearMesh( bool keepMemory = true );

	void createVertices( size_t count )
	{
		vertices.resize( count );
	}
	void createTriangles( size_t count )
	{
		triangles.resize( count );
	}

	// Get constant pointers to all vertices of the specified triangle
	void getTriangleVertices( uint32_t tri, const GEO::vec3** p1, const GEO::vec3** p2, const GEO::vec3** p3 ) const
	{
		const GEO::vec3i& t = getTriangle( tri );
		*p1 = &getVertex( t.x );
		*p2 = &getVertex( t.y );
		*p3 = &getVertex( t.z );
	}

	const GEO::vec3& getFirstTriangleVertex( uint32_t tri ) const
	{
		const GEO::vec3i& t = getTriangle( tri );
		return getVertex( t.x );
	}

	template<class Lambda>
	void generateVertices( uint32_t count, Lambda lambda )
	{
		vertices.resize( count );
		GEO::vec3* rdi = vertices.data();
		for( uint32_t i = 0; i < count; i++, rdi++ )
			*rdi = lambda( i );
	}

	template<class Lambda>
	void generateTriangles( uint32_t count, Lambda lambda )
	{
		triangles.resize( count );
		GEO::vec3i* rdi = triangles.data();
		for( uint32_t i = 0; i < count; i++, rdi++ )
			*rdi = lambda( i );
	}
};