#pragma once
#include "MeshBase.h"

// Tetrahedral volumetric mesh in 3D
class TetrahedralMesh : public MeshBase
{
	std::vector<__m128i> elements;

  public:
	void assignVertices( size_t count, const double* vb );

	// Set index buffer of the mesh
	void assignElements( size_t count, const uint32_t* ib );

	size_t countElements() const
	{
		return elements.size();
	}

	__m128i getElement( size_t elt ) const
	{
		return elements[ elt ];
	}

	void reorderMorton();

	// Compute 4.0 * specified coordinate of the center of the element; internal method for reorderMorton() implementation
	template<int COORD>
	inline double elementCenterX4( uint32_t idxElement ) const;
};