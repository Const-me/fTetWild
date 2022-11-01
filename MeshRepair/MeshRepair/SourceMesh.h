#pragma once
#include "API/library.h"
#include <floattetwild/Mesh.hpp>

namespace MeshRepair
{
	class SourceMesh : public ComLight::ObjectRoot<iSourceMesh>
	{
	public:
		HRESULT createMesh( uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib );

		floatTetWild::Mesh mesh;
	};
}  // namespace MeshRepair