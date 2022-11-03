#pragma once
#include "MeshBase.h"

// Tetrahedral volumetric mesh in 3D
class TetrahedralMesh : public MeshBase
{
	std::vector<__m128i> elements;

  public:

	// Compute 4.0 * specified coordinate of the center of the triangle
	template<int COORD>
	inline double elementCenterX4( uint32_t idxCell ) const;

	size_t countElements() const
	{
		return elements.size();
	}

	void reorderMorton();
};