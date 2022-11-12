#include "stdafx.h"
#include "VertexConnectedTets.h"
#include "LocalOperations.h"

#define USE_ORIGINAL_VERSION 0
#define COMPARE_VERSION 0

namespace
{
	struct alignas( 64 ) IntersectionBuffers
	{
		std::vector<int> a, b, c;
	};
	thread_local IntersectionBuffers buffers;
}  // namespace

#if CONN_TETS_SORTED_COPY
const std::vector<int>& floatTetWild::VertexConnectedTets::makeSortedVector() const
{
	if( !vec.empty() )
	{
		if( !sortedVec.empty() )
			return sortedVec;

		sortedVec = vec;
		std::sort( sortedVec.begin(), sortedVec.end() );
		return sortedVec;
	}

	sortedVec.clear();
	return sortedVec;
}
#endif

void floatTetWild::setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result )
{
#if USE_ORIGINAL_VERSION
	set_intersection( a.vec, b.vec, result );
#else
#if CONN_TETS_SORTED_COPY
	const std::vector<int>& as = a.makeSortedVector();
	const std::vector<int>& bs = b.makeSortedVector();
	std::set_intersection( as.begin(), as.end(), bs.begin(), bs.end(), std::back_inserter( result ) );
#else
	IntersectionBuffers& ib = buffers;
	ib.a = a.vec;
	ib.b = b.vec;
	std::sort( ib.a.begin(), ib.a.end() );
	std::sort( ib.b.begin(), ib.b.end() );
	std::set_intersection( ib.a.begin(), ib.a.end(), ib.b.begin(), ib.b.end(), std::back_inserter( result ) );
#endif
#endif
}

namespace
{
	inline void merge3_scalar( const int* a1, const int* a2, const int* b1, const int* b2, const int* c1, const int* c2, std::vector<int>& result )
	{
		int a = *a1, b = *b1, c = *c1;
		while( true )
		{
			if( a < b || a < c )
			{
				a1++;
				if( a1 >= a2 )
					return;
				a = *a1;
				continue;
			}
			if( b < a || b < c )
			{
				b1++;
				if( b1 >= b2 )
					return;
				b = *b1;
				continue;
			}
			if( c < a || c < b )
			{
				c1++;
				if( c1 >= c2 )
					return;
				c = *c1;
				continue;
			}
			result.push_back( a );
			a1++;
			b1++;
			c1++;
			if( a1 >= a2 || b1 >= b2 || c1 >= c2 )
				return;
			a = *a1;
			b = *b1;
			c = *c1;
		}
	}

	__forceinline __m128i load3( const int* x, const int* y, const int* z )
	{
		__m128i v = _mm_cvtsi32_si128( *x );
		v = _mm_insert_epi32( v, *y, 1 );
		v = _mm_insert_epi32( v, *z, 2 );
		return v;
	}

	inline void merge3_sse2( const int* a1, const int* a2, const int* b1, const int* b2, const int* c1, const int* c2, std::vector<int>& result )
	{
		while( true )
		{
			__m128i vec = load3( a1, b1, c1 );
			while( true )
			{
				const __m128i perm_baa = _mm_shuffle_epi32( vec, _MM_SHUFFLE( 3, 0, 0, 1 ) );
				const __m128i perm_ccb = _mm_shuffle_epi32( vec, _MM_SHUFFLE( 3, 1, 2, 2 ) );

				const __m128i lt1 = _mm_cmplt_epi32( vec, perm_baa );	 // Compare for [ a, b, c ] < [ b, a, a ]
				const __m128i lt2 = _mm_cmplt_epi32( vec, perm_ccb );	 // Compare for [ a, b, c ] < [ c, c, b ]
				const __m128i lt = _mm_or_si128( lt1, lt2 );			 // Combine with bitwise OR
				const uint32_t bmp = (uint32_t)_mm_movemask_epi8( lt );	 // Move bitmap to scalar register
				if( 0 != ( bmp & 1 ) )
				{
					// a < b || a < c
					a1++;
					if( a1 >= a2 )
						return;
					vec = _mm_insert_epi32( vec, *a1, 0 );
					continue;
				}
				if( 0 != ( bmp & 0x10 ) )
				{
					// b < a || b < c
					b1++;
					if( b1 >= b2 )
						return;
					vec = _mm_insert_epi32( vec, *b1, 1 );
					continue;
				}
				if( 0 != ( bmp & 0x100 ) )
				{
					// c < a || c < b
					c1++;
					if( c1 >= c2 )
						return;
					vec = _mm_insert_epi32( vec, *c1, 2 );
					continue;
				}

				// None of the above, we have the same integer in all 3 lanes of the vector
				assert( _mm_cvtsi128_si32( vec ) == _mm_extract_epi32( vec, 1 ) );
				assert( _mm_cvtsi128_si32( vec ) == _mm_extract_epi32( vec, 2 ) );

				result.push_back( _mm_cvtsi128_si32( vec ) );
				a1++;
				b1++;
				c1++;
				if( a1 >= a2 || b1 >= b2 || c1 >= c2 )
					return;
				break;	// Break into the outer loop to reload the complete vector from these 3 source pointers
			}
		}
	}

	inline void merge3( const std::vector<int>& a, const std::vector<int>& b, const std::vector<int>& c, std::vector<int>& result )
	{
		// merge3_scalar( a.data(), a.data() + a.size(), b.data(), b.data() + b.size(), c.data(), c.data() + c.size(), result );
		merge3_sse2( a.data(), a.data() + a.size(), b.data(), b.data() + b.size(), c.data(), c.data() + c.size(), result );
	}
}  // namespace

void floatTetWild::setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result )
{
#if USE_ORIGINAL_VERSION
	set_intersection( a.vec, b.vec, c.vec, result );
#else
	if( a.empty() || b.empty() || c.empty() )
		return;

#if CONN_TETS_SORTED_COPY
	const std::vector<int>& as = a.makeSortedVector();
	const std::vector<int>& bs = b.makeSortedVector();
	const std::vector<int>& cs = c.makeSortedVector();
	merge3( as, bs, cs, result );
#else
	IntersectionBuffers& ib = buffers;
	ib.a = a.vec;
	ib.b = b.vec;
	ib.c = c.vec;
	std::sort( ib.a.begin(), ib.a.end() );
	std::sort( ib.b.begin(), ib.b.end() );
	std::sort( ib.c.begin(), ib.c.end() );

	// std::set_intersection( ib.a.begin(), ib.a.end(), ib.b.begin(), ib.b.end(), std::back_inserter( result ) );
	// auto it = std::set_intersection( result.begin(), result.end(), ib.c.begin(), ib.c.end(), result.begin() );
	// result.resize( it - result.begin() );
	merge3( ib.a, ib.b, ib.c, result );
#endif

#if COMPARE_VERSION
	std::vector<int> orig;
	set_intersection( a.vec, b.vec, c.vec, orig );
	if( result != orig )
		__debugbreak();
#endif
#endif
}