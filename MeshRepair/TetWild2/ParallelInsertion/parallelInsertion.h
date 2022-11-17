#pragma once
#include "../src/Types.hpp"
#include <functional>

namespace floatTetWild
{
	using pfnInsertTri = std::function<void( int )>;

	void parallelInsertion(
	  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, const pfnInsertTri& pfn );
}  // namespace floatTetWild