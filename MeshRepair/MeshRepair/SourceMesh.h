#pragma once
#include "API/library.h"
#include "../TetWild2/src/Mesh.hpp"
#include "../TetWild2/Utils/Geogram2.h"

namespace MeshRepair
{
	class SourceMesh : public ComLight::ObjectRoot<iSourceMesh>
	{
	public:
		HRESULT createMesh( uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib );

		GEO2::Mesh mesh;

		// Another copy of the same data, stored in slightly different way.
		// TODO [RAM usage]: refactor away somehow
		std::vector<floatTetWild::Vector3> input_vertices;
		std::vector<floatTetWild::Vector3i> input_faces;
		std::vector<int> input_tags;
	};
}  // namespace MeshRepair