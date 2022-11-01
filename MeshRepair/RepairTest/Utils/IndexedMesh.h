#pragma once

// Indexed triangle mesh in system memory
struct IndexedMesh
{
	std::vector<DirectX::XMFLOAT3> vertices;
	std::vector<std::array<uint32_t, 3>> triangles;

	// Load binary STL from file on disk
	HRESULT loadBinaryStl( LPCTSTR path );

	// Save this mesh as a binary STL
	HRESULT saveBinaryStl( LPCTSTR path ) const;
};