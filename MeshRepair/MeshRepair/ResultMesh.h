#pragma once
#include "API/interfaces.cl.h"
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

	class ResultMesh32 : public ComLight::ObjectRoot<iResultMesh>
	{
		HRESULT ensureGoodMesh() const;

		HRESULT COMLIGHTCALL getSize( uint32_t& vertices, uint32_t& triangles ) const noexcept override final;
		HRESULT COMLIGHTCALL getIndexBuffer( uint32_t* rdi ) const noexcept override final;
		HRESULT COMLIGHTCALL getVertexBufferFP32( float* rdi ) const noexcept override final;
		HRESULT COMLIGHTCALL getVertexBufferFP64( double* rdi ) const noexcept override final;

	  public:
		std::vector<std::array<float, 3>> vertexBuffer;
		std::vector<std::array<uint32_t, 3>> indexBuffer;
	};
}  // namespace MeshRepair