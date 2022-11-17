#pragma once
#include "../src/Types.hpp"

namespace floatTetWild
{
	// Utility class to split set of triangles into disjoints sets, for parallel insertion
	class TrianglePartition
	{
		// True to preserve original order of the triangles
		// False to screw up the order, and treat output vectors as unordered sets. This last option is faster, avoids a few std::sort function calls.
		static constexpr bool preserveOrder = true;

		bool tryPartitionOnCoordinate(
		  const double* vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, double clearance, std::array<std::vector<int>, 3>& result );

		struct PartitionEntry
		{
			double min, max;
			uint32_t idxFaceInSet;
			uint32_t idxFaceInMesh;

			inline void computeMinMax( const double* vb, const std::vector<Vector3i>& ib, int triIndex );
		};
		std::vector<PartitionEntry> entries;
		struct SortEntries;

		void buildEntries( const double* vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces );

		bool tryPartitionEntries( size_t& leftEnd, size_t& middleEnd, double clearance ) const;

		struct EntryTemp
		{
			uint32_t idxFaceInSet;
			uint32_t idxFaceInMesh;
		};
		std::vector<EntryTemp> entriesTemp;
		struct SortTemp;

		void makePartitions( size_t leftEnd, size_t middleEnd, std::array<std::vector<int>, 3>& result );
		void makePartition( size_t begin, size_t end, std::vector<int>& dest );

	  public:
		// Try to partition the set of triangles into 3 subset, so the first 2 of these subsets can be processed in parallel
		// If successful, the first two sets of triangles are guaranteed to be separated by at least `clearance` distance.
		// If failed, the method returns false and does not change the output arrays
		bool tryPartition( const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance,
		  std::array<std::vector<int>, 3>& result );
	};
}  // namespace floatTetWild