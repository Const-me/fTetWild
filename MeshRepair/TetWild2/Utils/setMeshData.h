#pragma once
#include "../src/Types.hpp"
#include "Geogram2.h"

// Set vertex buffer of the mesh, upcasting coordinates to FP64
HRESULT assignMeshVertices( GEO2::Mesh& mesh, size_t count, const float* vb );

// Set index buffer of the mesh
HRESULT assignMeshTriangles( GEO2::Mesh& mesh, size_t count, const uint32_t* ib );

// Copy vertex and index buffer out of the mesh
HRESULT copyMeshData( const GEO2::Mesh& mesh, std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib );