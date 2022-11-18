#pragma once
#include "../src/Types.hpp"
#include <functional>

namespace floatTetWild
{
	using pfnInsertTri = std::function<void( int )>;

	// This version only uses 2 threads
	void parallelInsertionTwo(
	  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, const pfnInsertTri& pfn );

	inline void parallelInsertion(
	  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, const pfnInsertTri& pfn )
	{
		parallelInsertionTwo( vb, ib, faces, clearance, pfn );
	}
}  // namespace floatTetWild