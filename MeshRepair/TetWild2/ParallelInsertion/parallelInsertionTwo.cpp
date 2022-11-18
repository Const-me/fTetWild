#include "parallelInsertion.h"
#include "TrianglePartition.h"

namespace
{
	using namespace floatTetWild;

	static void insertSequential( const std::vector<int>& faces, const pfnInsertTri& pfn )
	{
		for( int i : faces )
			pfn( i );
	}

	static void insertOmp( const std::array<std::vector<int>, 3>& arr, const pfnInsertTri& pfn )
	{
#pragma omp parallel for schedule( dynamic, 1 )
		for( int s = 0; s < 2; s++ )
		{
			for( int i : arr[ s ] )
				pfn( i );
		}
		for( int i : arr[ 2 ] )
			pfn( i );
	}
}  // namespace

namespace floatTetWild
{
	void parallelInsertionTwo(
	  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, const pfnInsertTri& pfn )
	{
		// TODO: implement some recursion here, to use more than just 2 threads
		TrianglePartition part;
		std::array<std::vector<int>, 3> arr;
		if( part.tryPartition( vb, ib, faces, clearance, arr ) )
			insertOmp( arr, pfn );
		else
			insertSequential( faces, pfn );
	}
}  // namespace floatTetWild