#include <stdafx.h>
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
			const int* rsi = (const int*)buffer;
			const size_t columnHeight = len * 4;
			size_t columnStart = i0 * 4;
			const size_t columnsEnd = columnStart + columnSize * 4;
			for( ; columnStart < columnsEnd; columnStart += columnSize, rsi += columnHeight )
			{
				if( columnStart >= outputLength )
					continue;
				size_t columnEnd = columnStart + columnHeight;
				// The output buffer can be smaller than the input due to the INT_MAX values added to make the temporary buffer a multiple of 16 elements
				// That's why clipping the slices
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

	static __declspec( noinline ) __m128i* makeTempBuffer( size_t countBlocks )
	{
		std::vector<__m128i>& tmp = temporaryBuffer;
		tmp.resize( countBlocks * 4 );
		return tmp.data();
	}

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

			if( 0 != ( inputIntegers % 16 ) )
			{
				// Horizontally sort remaining 0-3 complete vectors
				const size_t completeVectors = ( inputIntegers % 16 ) / 4;
				__m128i* const rdiEnd = rdi + 4;
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
				for( ; rdi < rdiEnd; rdi++ )
					_mm_storeu_si128( rdi, Which::maxVector() );
			}

			return countBlocks * 4;
		}

		static __forceinline size_t divByShrinkFactor( size_t n )
		{
			// 5.1. Implementation details "The shrink factor for our in core sorting algorithm was 1.28"
			// Probably not coincidentally, the inverse of that value is 25/32, the product can be computed with two fast instructions
			return ( n * 25 ) / 32;
		}

		// Implement combsort on the array into the transposed order
		static void sortStep2( __m128i* buffer, size_t countVectors )
		{
			__m128i* const ptr = buffer;

			for( size_t gap = divByShrinkFactor( countVectors ); gap > 1; gap = divByShrinkFactor( gap ) )
			{
				// Straight comparisons
				for( size_t i = 0; i < countVectors - gap; i++ )
					Which::vector_cmpswap( ptr + i, ptr + i + gap );

				// Skewed comparisons, when i + gap exceeds N/4
				for( size_t i = countVectors - gap; i < countVectors; i++ )
					Which::vector_cmpswap_skew( ptr + i, ptr + i + gap - countVectors );
			}

			do
			{
				for( size_t i = 0; i < countVectors - 1; i++ )
					Which::vector_cmpswap( ptr + i, ptr + i + 1 );
				Which::vector_cmpswap_skew( ptr + countVectors - 1, ptr );
			} while( !Which::isSorted( ptr, countVectors ) );
		}

	  public:
		static void sort( const std::vector<E>& source, std::vector<E>& dest )
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
				buffer = makeTempBuffer( countBlocks );
				dest.resize( source.size() );
			}

			// ==== Step #1 of the algorithm ====
			// move to the transposed vector, while sorting lanes in the vector
			// Also pad the remainder with INT_MAX / UINT_MAX
			const size_t countVectors = sortStep1( source, buffer, countBlocks );

			// ==== Step #2 of the algorithm, vertically sort these vectors ====
			sortStep2( buffer, countVectors );

			// ==== Step #3 of the algorithm, transpose the output back to normal ====
			if( inPlace )
			{
				transposeInplace( buffer, countVectors / 4 );
				// Truncate vector to the original size
				dest.resize( source.size() );
			}
			else
				transposeWithCopy( buffer, countVectors / 4, (int*)dest.data(), dest.size() );
		}

		static void sortInPlace( std::vector<E>& vec )
		{
			const size_t sourceSize = vec.size();
			if( 0 == sourceSize )
				return;

			const size_t countBlocks = ( sourceSize + 15 ) / 16;
			const bool inPlace = countBlocks <= maxInplaceBlocks();
			__m128i* buffer;
			if( inPlace )
			{
				// Slightly expand the vector to the next multiple of 16 elements, appending maximum value for the integer type
				vec.resize( countBlocks * 16, Which::maxValue );
				buffer = (__m128i*)vec.data();

				// Note the first step of the algorithm becomes way more efficient
				__m128i* const bufferEnd = buffer + countBlocks * 4;
				for( __m128i* p = buffer; p < bufferEnd; p += 4 )
					Which::copySortedBlock( p, p );
			}
			else
			{
				// Unfortunately, we need a local copy for the large final transpose
				// We only have a lookup table for up to 64 4x4 blocks in the vector
				buffer = makeTempBuffer( countBlocks );
				sortStep1( vec, buffer, countBlocks );
			}

			// ==== Step #2 of the algorithm, execute combsort on the vector integer array into the transposed order ====
			sortStep2( buffer, countBlocks * 4 );

			// ==== Step #3 of the algorithm, transpose the output back to normal ====
			if( inPlace )
			{
				transposeInplace( buffer, countBlocks );
				// Truncate vector to the original size
				vec.resize( sourceSize );
			}
			else
				transposeWithCopy( buffer, countBlocks, (int*)vec.data(), sourceSize );
		}
	};
}  // namespace

namespace AASort
{
	void sortVector( const std::vector<int>& source, std::vector<int>& dest )
	{
		SortImpl<int, Primitives<IntTraits>>::sort( source, dest );
	}
	void sortVector( std::vector<int>& vec )
	{
		SortImpl<int, Primitives<IntTraits>>::sortInPlace( vec );
	}

	void sortVector( const std::vector<uint32_t>& source, std::vector<uint32_t>& dest )
	{
		SortImpl<uint32_t, Primitives<UintTraits>>::sort( source, dest );
	}
	void sortVector( std::vector<uint32_t>& vec )
	{
		SortImpl<uint32_t, Primitives<UintTraits>>::sortInPlace( vec );
	}
};	// namespace AASort