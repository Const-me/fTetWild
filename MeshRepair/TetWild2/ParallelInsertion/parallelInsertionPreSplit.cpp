#include "parallelInsertion.h"
#include "TrianglePartition.h"

namespace
{
	using namespace floatTetWild;

	class PartitionsPool
	{
		std::vector<std::vector<int>> freedElemens;

	  public:
		std::vector<int> emplace()
		{
			if( !freedElemens.empty() )
			{
				std::vector<int> res = std::move( freedElemens.back() );
				freedElemens.pop_back();
				res.clear();
				return res;
			}
			else
			{
				return std::vector<int> {};
			}
		}

		void release( std::vector<int>& vec )
		{
			freedElemens.emplace_back( std::move( vec ) );
		}
	};

	using NestedVectors = std::vector<std::vector<std::vector<int>>>;

	struct SplitContext
	{
		const std::vector<Vector3>& vb;
		const std::vector<Vector3i>& ib;
		const __m256d clearance;
		NestedVectors& vectors;

		PartitionsPool pool;
		TrianglePartition part;

		void recursion( std::vector<int>& faces, size_t level, bool root )
		{
			std::array<std::vector<int>, 3> arr { pool.emplace(), pool.emplace(), pool.emplace() };

			if( part.tryPartition( vb, ib, faces, clearance, arr ) )
			{
				if( !root )
				{
					// We no longer need the input array; unless it's the parameter to parallelInsertionPreSplit(), release to the pool
					pool.release( faces );
				}
				recursion( arr[ 0 ], level + 1, false );
				recursion( arr[ 1 ], level + 1, false );
				if( !arr[ 2 ].empty() )
				{
					// We can't compute different pieces of the arr[ 2 ] in parallel with arr[ 0 ] / arr[ 1 ], they gonna modify same elements of the mesh
					// It is possible to implement this better, though.
					// Need to build a tree of these dependent pieces, then extract parallelizeable chunks from the tree
					addToLevel( level, arr[ 2 ], false );
				}
				else
					pool.release( arr[ 2 ] );
			}
			else
			{
				// Partition failed, release all 3 temporary arrays to the pool
				for( ptrdiff_t i = 2; i >= 0; i-- )
					pool.release( arr[ i ] );
				// Schedule the whole array at the current level
				addToLevel( level, faces, root );
			}
		}

		void addToLevel( size_t level, std::vector<int>& faces, bool root )
		{
			if( level >= vectors.size() )
				vectors.resize( level + 1 );
			if( root )
				vectors[ level ].push_back( faces );
			else
				vectors[ level ].emplace_back( std::move( faces ) );
		}
	};

	static void preSplit(
	  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, NestedVectors& vectors )
	{
		SplitContext context { vb, ib, clearance, vectors };
		context.recursion( const_cast<std::vector<int>&>( faces ), 0, true );
	}

	static void insertionOmp( const NestedVectors& vectors, const pfnInsertTri& pfn )
	{
		for( size_t lvl = vectors.size() - 1; lvl > 0; lvl-- )
		{
			const auto& level = vectors[ lvl ];
			const int len = (int)level.size();
#pragma omp parallel for schedule( dynamic, 1 )
			for( int s = 0; s < len; s++ )
			{
				for( int t : level[ s ] )
					pfn( t );
			}
		}
		if( !vectors[ 0 ].empty() )
			for( int i : vectors[ 0 ][ 0 ] )
				pfn( i );
	}
}  // namespace

namespace floatTetWild
{
	void __declspec( noinline ) parallelInsertionPreSplit(
	  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, const pfnInsertTri& pfn )
	{
		NestedVectors vectors;
		preSplit( vb, ib, faces, clearance, vectors );
		assert( vectors[ 0 ].size() <= 1 );

		if( vectors.size() > 1 )
			insertionOmp( vectors, pfn );
		else
		{
			assert( vectors.size() == 1 );
			for( int i : vectors[ 0 ][ 0 ] )
				pfn( i );
		}
	}
}  // namespace floatTetWild