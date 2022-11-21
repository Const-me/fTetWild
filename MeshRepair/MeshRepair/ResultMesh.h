#pragma once
#include "API/library.h"
#include "../ComLightLib/comLightServer.h"
#include "../TetWild2/src/Types.hpp"

namespace MeshRepair
{
	class ResultMesh : public ComLight::ObjectRoot<iResultMesh>
	{
		HRESULT ensureGoodMesh() const;

		HRESULT COMLIGHTCALL getSize( uint32_t& vertices, uint32_t& triangles ) const noexcept override final;
		HRESULT COMLIGHTCALL getIndexBuffer( uint32_t* rdi ) const noexcept override final;
		HRESULT COMLIGHTCALL getVertexBufferFP32( float* rdi ) const noexcept override final;
		HRESULT COMLIGHTCALL getVertexBufferFP64( double* rdi ) const noexcept override final;

	  public:
		floatTetWild::VertexBuffer V_sf;
		floatTetWild::SurfaceIndexBuffer F_sf;
	};
}  // namespace MeshRepair