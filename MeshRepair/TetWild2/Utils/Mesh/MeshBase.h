#pragma once
#include "../../src/Types.hpp"
#include "../Geogram2.h"

// Base class for meshes, contains just the vertex buffer, with 3D coordinates in FP64 precision
class MeshBase
{
  protected:
	std::vector<GEO2::vec3> vertices;

  public:
	  using vec3 = GEO2::vec3;

	// Set vertex buffer of the mesh, upcasting coordinates to FP64
	HRESULT assignVertices( size_t count, const float* vb );

	uint32_t countVertices() const
	{
		return vertices.size();
	}

	const vec3* vertexPointer() const
	{
		return vertices.data();
	}
	vec3* vertexPointer()
	{
		return vertices.data();
	}
	const vec3& getVertex( uint32_t v ) const
	{
		return vertices[ v ];
	}
	void createVertices( size_t count )
	{
		vertices.resize( count );
	}

	template<class Lambda>
	void generateVertices( uint32_t count, Lambda lambda )
	{
		vertices.resize( count );
		vec3* rdi = vertices.data();
		for( uint32_t i = 0; i < count; i++, rdi++ )
			*rdi = lambda( i );
	}
};