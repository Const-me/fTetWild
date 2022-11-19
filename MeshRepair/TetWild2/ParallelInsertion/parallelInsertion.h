#pragma once
#include "../src/Types.hpp"
#include <functional>

namespace floatTetWild
{
	using pfnInsertTri = std::function<void( int )>;

	// The initial version which only uses 2 threads
	void parallelInsertionTwo(
	  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, const pfnInsertTri& pfn );

	// This version pre-splits triangles on the main thread, then runs the main loops using OpenMP
	void parallelInsertionPreSplit(
	  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, const pfnInsertTri& pfn );

	// Run insertion in parallel, when possible
	inline void parallelInsertion(
	  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, const pfnInsertTri& pfn )
	{
		// parallelInsertionTwo( vb, ib, faces, clearance, pfn );
		parallelInsertionPreSplit( vb, ib, faces, clearance, pfn );
	}
}  // namespace floatTetWild