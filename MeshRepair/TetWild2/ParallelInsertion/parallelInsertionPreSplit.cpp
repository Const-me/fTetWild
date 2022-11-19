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
					pool.release( faces );
				recursion( arr[ 0 ], level + 1, false );
				recursion( arr[ 1 ], level + 1, false );
				recursion( arr[ 2 ], level, false );
			}
			else
			{
				for( ptrdiff_t i = 2; i >= 0; i-- )
					pool.release( arr[ i ] );
				addToLevel( level, faces );
			}
		}

		void addToLevel( size_t level, std::vector<int>& faces )
		{
			if( level >= vectors.size() )
				vectors.resize( level + 1 );
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
		assert( vectors[ 0 ].size() == 1 );

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