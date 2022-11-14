#pragma once
#include <array>
#include <immintrin.h>

namespace floatTetWild
{
	struct alignas( 16 ) ElementInQueue
	{
		std::array<int, 2> v_ids;
		double weight;

		ElementInQueue() = default;
		ElementInQueue( const std::array<int, 2>& ids, double w )
			: v_ids( ids )
			, weight( w )
		{
		}
		ElementInQueue( int e0, int e1, double w )
			: v_ids { e0, e1 }
			, weight( w )
		{
		}

		static inline bool compareElements( const ElementInQueue& e1, const ElementInQueue& e2 )
		{
			// The code looks scary, but this function only compiles to 11 instructions
			// Even better, most of them run in parallel, there're 2 parallel dependency chains in that code.
			// Half of the intrinsics below are casting types to workaround C++ type safety, they're free in runtime
			static_assert( sizeof( ElementInQueue ) == 16 );

			// Load both values to vector registers
			__m128d v1 = _mm_load_pd( (const double*)&e1 );
			__m128d v2 = _mm_load_pd( (const double*)&e2 );

			// Compare FP64 lanes for both e1 > e2 and e1 < e2
			__m128 gtw = _mm_castpd_ps( _mm_cmpgt_pd( v1, v2 ) );
			__m128 ltw = _mm_castpd_ps( _mm_cmplt_pd( v1, v2 ) );

			// Compare int32 lanes for both e1 > e2 and e1 < e2
			__m128 gti = _mm_castsi128_ps( _mm_cmpgt_epi32( _mm_castpd_si128( v1 ), _mm_castpd_si128( v2 ) ) );
			__m128 lti = _mm_castsi128_ps( _mm_cmplt_epi32( _mm_castpd_si128( v1 ), _mm_castpd_si128( v2 ) ) );

			// The original code uses operator < of the std::array<int, 2> class
			// Apparently, that class treats first integer as the most significant result, similar to the way memcmp treats bytes
			// To compensate, flipping XY lanes while combining these vectors. Otherwise, we could do that with _mm_blend_ps, slightly faster
			// Also note we merge different comparisons - because when weights are equal, the comparison of the integers is inverted
			constexpr int shuff = _MM_SHUFFLE( 3, 2, 0, 1 );
			__m128 gt = _mm_shuffle_ps( lti, gtw, shuff );
			__m128 lt = _mm_shuffle_ps( gti, ltw, shuff );

			uint32_t maskGt = (uint32_t)_mm_movemask_ps( gt );
			uint32_t maskLt = (uint32_t)_mm_movemask_ps( lt );
			return maskGt > maskLt;
		}
	};

	struct cmp_l
	{
		bool compareOld( const ElementInQueue& e1, const ElementInQueue& e2 )
		{
			if( e1.weight == e2.weight )
				return e1.v_ids > e2.v_ids;
			return e1.weight < e2.weight;
		}

		bool operator()( const ElementInQueue& e1, const ElementInQueue& e2 )
		{
#if 1
			return ElementInQueue::compareElements( e2, e1 );
#else
			bool scalar = compareOld( e1, e2 );
			bool vec = ElementInQueue::compareElements( e2, e1 );
			if( scalar != vec )
				__debugbreak();
			return vec;
#endif
		}
	};

	struct cmp_s
	{
		bool compareOld( const ElementInQueue& e1, const ElementInQueue& e2 )
		{
			if( e1.weight == e2.weight )
				return e1.v_ids < e2.v_ids;
			return e1.weight > e2.weight;
		}

		bool operator()( const ElementInQueue& e1, const ElementInQueue& e2 )
		{
#if 1
			return ElementInQueue::compareElements( e1, e2 );
#else
			bool scalar = compareOld( e1, e2 );
			bool vec = ElementInQueue::compareElements( e1, e2 );
			if( scalar != vec )
				__debugbreak();
			return vec;
#endif
		}
	};
}  // namespace floatTetWild