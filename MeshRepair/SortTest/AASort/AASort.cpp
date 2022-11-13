#include "AASort.h"
#include "AASortUtils.hpp"
#include "AASortTranspose.h"

namespace
{
	using namespace AASort;

	// First pass of the transpose, process 4x4 matrices
	__forceinline void transpose4x4Blocks( __m128i* const begin, size_t countBlocks )
	{
		__m128i* const rdiEnd = begin + countBlocks * 4;

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
	}

	// Step #3 of the algorithm, when vector size does not exceed 4kb = 1024 elements
	static void transposeInplace( __m128i* const begin, size_t countBlocks )
	{
		transpose4x4Blocks( begin, countBlocks );
		transposeBlocksInplace( begin, countBlocks );
	}

	static void transposeWithCopy( __m128i* const begin, const size_t countBlocks, int* const rdi, const size_t outputLength )
	{
		assert( countBlocks == ( ( outputLength + 15 ) / 16 ) );
		const size_t blockLength = maxInplaceBlocks();
		const size_t outerLoopLength = ( countBlocks + blockLength - 1 ) / blockLength;

		// Integers per column in the complete buffer
		const size_t columnSize = countBlocks * 4;

		for( size_t i = 0; i < outerLoopLength; i++ )
		{
			const size_t i0 = i * blockLength;
			const size_t i2 = std::min( ( i + 1 ) * blockLength, countBlocks );
			assert( i2 > i0 );

			// Transpose a middle portion of the buffer in-place. That portion is up to 4kb, fits in L1D cache
			__m128i* const buffer = begin + i0 * 4;
			const size_t len = i2 - i0;
			transposeInplace( buffer, len );

			// Copy 4 slices of the newly transposed portion into the corresponding location of the output buffer
			// The output buffer can be smaller than the input due to the INT_MAX values added to make the temporary buffer a multiple of 16 elements
			// That's why we need to clip these slices
			const int* rsi = (const int*)buffer;
			const size_t columnHeight = len * 4;
			size_t columnStart = i0 * 4;
			const size_t columnsEnd = columnStart + columnSize * 4;
			for( ; columnStart < columnsEnd; columnStart += columnSize, rsi += columnHeight )
			{
				if( columnStart >= outputLength )
					continue;
				size_t columnEnd = columnStart + columnHeight;
				columnEnd = std::min( columnEnd, outputLength );
				if( columnEnd <= columnStart )
					continue;
				memcpy( rdi + columnStart, rsi, ( columnEnd - columnStart ) * 4 );
			}
		}
	}

	// This temporary buffer is only used when input vector exceeds 1024 elements
	// Otherwise, the algorithm works in-place
	static thread_local std::vector<__m128i> temporaryBuffer;

	template<class E, class Which>
	class SortImpl
	{
		static size_t sortStep1( const std::vector<E>& source, __m128i* rdi, size_t countBlocks )
		{
			const size_t inputIntegers = source.size();
			assert( countBlocks == ( ( inputIntegers + 15 ) / 16 ) );

			const E* rsi = source.data();
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

	  public:
		static __forceinline void sort( const std::vector<E>& source, std::vector<E>& dest )
		{
			if( source.empty() )
			{
				dest.clear();
				return;
			}

			const size_t countBlocks = ( source.size() + 15 ) / 16;
			const bool inPlace = countBlocks <= maxInplaceBlocks();
			__m128i* buffer;
			if( inPlace )
			{
				dest.resize( countBlocks * 16 );
				buffer = (__m128i*)dest.data();
			}
			else
			{
				temporaryBuffer.resize( countBlocks * 4 );
				dest.resize( source.size() );
				buffer = temporaryBuffer.data();
			}

			// ==== Step #1 of the algorithm ====
			// move to the transposed vector, while sorting lanes in the vector
			// Also pad the remainder with INT_MAX / UINT_MAX
			const size_t countVectors = sortStep1( source, buffer, countBlocks );

			// ==== Step #2 of the algorithm, vertically sort these vectors ====
			__m128i* const ptr = buffer;
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
			if( inPlace )
			{
				transposeInplace( buffer, countVectors / 4 );
				dest.resize( source.size() );
			}
			else
				transposeWithCopy( buffer, countVectors / 4, (int*)dest.data(), dest.size() );
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