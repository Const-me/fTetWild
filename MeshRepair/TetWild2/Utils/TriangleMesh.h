#pragma once
#include <geogram/basic/geometry.h>

class TriangleMesh
{
	GEO::vector<GEO::vec3> vertices;
	GEO::vector<GEO::vec3i> triangles;

  public:
	// returns 3.0 * specified coordinate of the center of the triangle
	template<int COORD>
	inline double triangleCenterX3( uint32_t idxTri ) const;

	// Set vertex buffer of the mesh, upcasting coordinates to FP64
	HRESULT assignVertices( size_t count, const float* vb );

	// Set index buffer of the mesh
	HRESULT assignTriangles( size_t count, const uint32_t* ib );

	size_t countTriangles() const
	{ 
		return triangles.size();
	}

	void reorderMorton();
};