#pragma once
#include "API/interfaces.cl.h"
#include "../ComLightLib/comLightServer.h"
#include "../TetWild2/src/Mesh.h"
#include "../TetWild2/Utils/Mesh/TriangleMesh.h"

namespace MeshRepair
{
	class SourceMesh : public ComLight::ObjectRoot<iSourceMesh>
	{
		GEO2::Mesh mesh;

	  public:

		HRESULT createMesh( uint32_t countVertices, const float* vb, uint32_t countTriangles, const uint32_t* ib );

		void makeBuffers( std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib, std::vector<int>* tags ) const;

		const GEO2::Mesh& getMesh() const
		{
			return mesh;
		}
	};
}  // namespace MeshRepair