#pragma once
#include "Types.hpp"

namespace floatTetWild
{
	constexpr int8_t NOT_SURFACE = 127;
	constexpr int8_t KNOWN_NOT_SURFACE = -63;
	constexpr int8_t KNOWN_SURFACE = 63;

	constexpr int8_t NO_SURFACE_TAG = 0;
	constexpr int8_t NOT_BBOX = -1;
	constexpr int OPP_T_ID_UNKNOWN = -2;
	constexpr int OPP_T_ID_BOUNDARY = -1;

	struct MeshTet
	{
		Vector4i indices;

		MeshTet()
		{
		}
		MeshTet( const Vector4i& idx )
			: indices( idx )
		{
		}
		MeshTet( int v0, int v1, int v2, int v3 )
			: indices( v0, v1, v2, v3 )
		{
		}

		inline void reset()
		{
			is_surface_fs = { { NOT_SURFACE, NOT_SURFACE, NOT_SURFACE, NOT_SURFACE } };
			is_bbox_fs = { { NOT_BBOX, NOT_BBOX, NOT_BBOX, NOT_BBOX } };
			opp_t_ids = { { OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN } };
			surface_tags = { { 0, 0, 0, 0 } };

			quality = 0;
			scalar = 0;
			is_removed = false;
			is_outside = false;
		}

		inline const Vector4i& pts_indices() const
		{
			return indices;
		}

		inline int& operator[]( const int index )
		{
			assert( index >= 0 && index < 4 );
			return indices[ index ];
		}

		inline int operator[]( const int index ) const
		{
			assert( index >= 0 && index < 4 );
			return indices[ index ];
		}

		// If the argument is in the indices vector, integer in [ 0 .. 3 ] interval with the position. Otherwise, -1
		inline int find( int ele ) const
		{
			static_assert( sizeof( indices ) == 16 );

			// Compare vector elements for ( e == ele )
			const __m128i vec = _mm_loadu_si128( (const __m128i*)indices.data() );
			const __m128i needle = _mm_set1_epi32( ele );
			const __m128i eq = _mm_cmpeq_epi32( needle, vec );

			// Move into bitmap in a scalar register, 0 = not equal, 1 = equal
			const uint32_t bmp = (uint32_t)_mm_movemask_ps( _mm_castsi128_ps( eq ) );

#ifdef __AVX__
			if( 0 != bmp )
				return (int)_tzcnt_u32( bmp );
			return -1;
#else
			unsigned long bitIndex;
			if( _BitScanForward( &bitIndex, bmp ) )
				return (int)bitIndex;
			return -1;
#endif
		}

		// Find first index not equal to either of the arguments
		inline int find_opp( int v0_id, int v1_id, int v2_id ) const
		{
			// Compare vector elements for ( e == v0 ) || ( e == v1 ) || ( e == v2 )
			const __m128i vec = _mm_loadu_si128( (const __m128i*)indices.data() );

			__m128i needle = _mm_set1_epi32( v0_id );
			__m128i eq = _mm_cmpeq_epi32( vec, needle );

			needle = _mm_set1_epi32( v1_id );
			eq = _mm_or_si128( eq, _mm_cmpeq_epi32( vec, needle ) );

			needle = _mm_set1_epi32( v2_id );
			eq = _mm_or_si128( eq, _mm_cmpeq_epi32( vec, needle ) );

			// Move into bitmap in a scalar register, 0 = not equal, 1 = equal
			uint32_t bmp = (uint32_t)_mm_movemask_ps( _mm_castsi128_ps( eq ) );
			// Invert lower 4 bits in the bitmap, results in 0 = equal to at least one of the inputs, 1 = not equal to either of the inputs
			bmp ^= 0b1111;

			// Scan for the lowest set bit in the bitmap
#ifdef __AVX__
			if( 0 != bmp )
				return (int)_tzcnt_u32( bmp );
			return -1;
#else
			unsigned long bitIndex;
			if( _BitScanForward( &bitIndex, bmp ) )
				return (int)bitIndex;
			return -1;
#endif
		}

		inline void print( const Logger& log, int id ) const
		{
			log.logDebug( "%i: [ %i, %i, %i, %i ]", id, indices[ 0 ], indices[ 1 ], indices[ 2 ], indices[ 3 ] );
		}

		std::array<int8_t, 4> is_surface_fs = { { NOT_SURFACE, NOT_SURFACE, NOT_SURFACE, NOT_SURFACE } };

		// For every face of the tetrahedron, this value tells whether the face is co-planar with some surface of the bounding box of the domain.
		// -1 means the face is inside, [ 0 .. 5 ] means the face is on the bounding box
		// The tetrahedron's faces are numbered [ 0, 1, 2 ], [ 1, 2, 3 ], [ 2, 3, 0 ], [ 3, 0, 1 ]
		// The bounding box faces are numbered -X, +X, -Y, +Y, -Z, +Z
		std::array<int8_t, 4> is_bbox_fs = { { NOT_BBOX, NOT_BBOX, NOT_BBOX, NOT_BBOX } };

		std::array<int8_t, 4> surface_tags = { { 0, 0, 0, 0 } };
		bool is_removed = false;
		bool is_outside = false;
		uint8_t scalar = 0;
		std::array<int, 4> opp_t_ids = { { OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN } };

		double quality = 0;
	};
}