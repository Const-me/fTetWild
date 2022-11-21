#pragma once
#include "../../TetWild2/src/Types.hpp"

void downcastMesh( const floatTetWild::VertexBuffer& vertices, const floatTetWild::SurfaceIndexBuffer& faces, std::vector<std::array<float, 3>>& vb32,
  std::vector<std::array<uint32_t, 3>>& ib32 );

void downcastFloats( float* rdi, size_t count, const double* rsi );