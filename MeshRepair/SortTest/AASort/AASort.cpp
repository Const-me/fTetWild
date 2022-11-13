#include "AASort.h"
#include "AASortUtils.hpp"
#include "AASortTranspose.h"

namespace
{
	static void transposeResult( __m128i* const begin, size_t countBlocks )
	{
		__m128i* const rdiEnd = begin + countBlocks * 4;

		// First pass, transpose 4x4 matrices
		for( __m128i* p = begin; p < rdiEnd; p += 4 )
		{
			__m128i v0 = _mm_loadu_si128( p );
			__m128i v1 = _mm_loadu_si128( p + 1 );
			__m128i v2 = _mm_loadu_si128( p + 2 );
			__m128i v3 = _mm_loadu_si128( p + 3 );

			transpose4x4( v0, v1, v2, v3 );

			_mm_storeu_si128( p, v0 );
			_mm_storeu_si128( p + 1, v1 );
			_mm_storeu_si128( p + 2, v2 );
			_mm_storeu_si128( p + 3, v3 );
		}

		// Second transpose pass, the implementation depends on the length of the vector
		// When the data size doesn't exceed 256 vectors = 4kb, TransposePermutations.inl file contains a lookup table for single-pass transposes
		// Larger vectors should do something else 
		AASort::transposeResultOuter( begin, countBlocks );
	}

	template<class E, class Which>
	class SortImpl
	{
		static __forceinline size_t sortStep1( const std::vector<E>& source, std::vector<E>& dest )
		{
			const size_t inputIntegers = source.size();
			if( 0 == inputIntegers )
			{
				dest.clear();
				return 0;
			}

			const size_t countBlocks = ( inputIntegers + 15 ) / 16;
			dest.resize( countBlocks * 16 );

			const E* rsi = source.data();
			__m128i* rdi = (__m128i*)dest.data();
			// Vertically sort and transpose complete 4x4 blocks
			const size_t completeBlocks = inputIntegers / 16;
			for( size_t i = 0; i < completeBlocks; i++ )
			{
				Which::copySortedBlock( rdi, rsi );
				rsi += 16;
				rdi += 4;
			}

			if( 0 == ( inputIntegers % 16 ) )
				return countBlocks * 4;

			// Horizontally sort remaining 0-3 complete vectors
			const size_t completeVectors = ( inputIntegers % 16 ) / 4;
			const __m128i* const rdiEnd = rdi + 4;
			for( size_t i = 0; i < completeVectors; i++ )
			{
				__m128i v = _mm_loadu_si128( (const __m128i*)rsi );
				rsi += 4;
				v = Which::sortLanes( v );
				_mm_storeu_si128( rdi, v );
				rdi++;
			}
			// Unless the input size is a multiple of 4, horizontally sort the incomplete vector
			const size_t incompleteVec = ( inputIntegers % 4 );
			if( 0 != incompleteVec )
			{
				__m128i v = Which::loadPartialMaxPad( rsi, incompleteVec );
				v = Which::sortLanes( v );
				_mm_storeu_si128( rdi, v );
				rdi++;
			}
			// Fill the remaining 0-3 destination vectors with maximum integer value
			for( ; rdi < rdiEnd; rsi++ )
			{
				__m128i v = Which::maxVector();
				_mm_storeu_si128( rdi, v );
				rdi++;
			}
			return countBlocks * 4;
		}

		static __forceinline size_t divByShrinkFactor( size_t n )
		{
			// The authors used 1.3 for the shrink factor
			// We don't want any FP math here, due to the latency
			return ( n * 10 ) / 13;
		}

		static __forceinline void transposeVector( std::vector<E>& dest )
		{
			assert( !dest.empty() );
			assert( 0 == ( dest.size() % 16 ) );

			const size_t countBlocks = dest.size() / 16;
			transposeResult( (__m128i*)dest.data(), countBlocks );
		}

	  public:
		static __forceinline void sort( const std::vector<E>& source, std::vector<E>& dest )
		{
			// ==== Step #1 of the algorithm ====
			// move to the transposed vector, while sorting lanes in the vector
			// Also pad the remainder with INT_MAX / UINT_MAX
			const size_t countVectors = sortStep1( source, dest );
			if( 0 == countVectors )
				return;

			// ==== Step #2 of the algorithm, vertically sort these vectors ====
			__m128i* const ptr = (__m128i*)dest.data();
			size_t gap = divByShrinkFactor( countVectors );

			while( gap > 1 )
			{
				// Straight comparisons
				for( size_t i = 0; i < countVectors - gap; i++ )
					Which::vector_cmpswap( ptr + i, ptr + i + gap );
				// Skewed comparisons, when i + gap exceeds N/4
				for( size_t i = countVectors - gap; i < countVectors; i++ )
					Which::vector_cmpswap_skew( ptr + i, ptr + i + gap - countVectors );

				// Shrink the gap
				gap = divByShrinkFactor( gap );
			}

			while( true )
			{
				for( size_t i = 0; i < countVectors - 1; i++ )
					Which::vector_cmpswap( ptr + i, ptr + i + 1 );

				Which::vector_cmpswap_skew( ptr + countVectors - 1, ptr );
				if( Which::isSorted( ptr, countVectors ) )
					break;
			}

			// ==== Step #3 of the algorithm, transpose the output back to normal ====
			transposeVector( dest );
			dest.resize( source.size() );
		}
	};
}  // namespace

namespace AASort
{
	void sortVector( const std::vector<int>& source, std::vector<int>& dest )
	{
		SortImpl<int, Primitives<IntTraits>>::sort( source, dest );
	}

	void sortVector( const std::vector<uint32_t>& source, std::vector<uint32_t>& dest )
	{
		SortImpl<uint32_t, Primitives<UintTraits>>::sort( source, dest );
	}
};	// namespace AASort