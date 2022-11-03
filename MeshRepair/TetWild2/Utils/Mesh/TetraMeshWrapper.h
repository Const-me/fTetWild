#pragma once
#include <geogram/mesh/mesh_AABB.h>

class TetraMeshWrapper : public GEO::Mesh
{
  public:
	GEO::Mesh& unwrap()
	{
		return *this;
	}
};

class MeshCellsAABBWrapper : public GEO::MeshCellsAABB
{
  public:
	MeshCellsAABBWrapper( TetraMeshWrapper& mesh, bool reorder = true )
		: GEO::MeshCellsAABB( mesh.unwrap(), reorder )
	{
	}
};