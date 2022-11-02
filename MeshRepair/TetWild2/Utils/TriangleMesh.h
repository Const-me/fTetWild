#pragma once
#include <geogram/basic/geometry.h>

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
};