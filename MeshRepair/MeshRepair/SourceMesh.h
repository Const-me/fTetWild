#pragma once
#include "API/library.h"
namespace MeshRepair
{
	class SourceMesh : public ComLight::ObjectRoot<iSourceMesh>
	{
	public:
		HRESULT createMesh( uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib );
	};
}  // namespace MeshRepair