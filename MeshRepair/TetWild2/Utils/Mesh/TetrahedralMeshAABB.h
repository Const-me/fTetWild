#pragma once
class TetrahedralMesh;

class TetrahedralMeshAABB
{
	struct alignas( 16 ) BBox
	{
		std::array<double, 3> min, max;
		inline bool containsPoint( const GEO::vec3& p ) const;
	};

	std::vector<BBox> boxes;
	const TetrahedralMesh& mesh;

	// Private methods for initialization
	void initializeRecursive( size_t node_index, size_t b, size_t e );
	void computeBox( BBox& rdi, size_t index ) const;
	static inline void mergeBoxes( BBox& rdi, const BBox& a, const BBox& b );

	uint32_t containingElement( const GEO::vec3& p, size_t n, size_t b, size_t e ) const;

  public:

	TetrahedralMeshAABB( TetrahedralMesh& sourceMesh, bool reorder = true );

	static constexpr uint32_t NO_TET = UINT_MAX;

	uint32_t containingElement( const GEO::vec3& p ) const;
};