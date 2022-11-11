#pragma once
#include "Types.hpp"

namespace floatTetWild
{
	constexpr int8_t NOT_SURFACE = 127;
	constexpr int8_t KNOWN_NOT_SURFACE = -63;
	constexpr int8_t KNOWN_SURFACE = 63;

#define NO_SURFACE_TAG 0
#define NOT_BBOX -1
#define OPP_T_ID_UNKNOWN -2
#define OPP_T_ID_BOUNDARY -1

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

		inline int find( int ele ) const
		{
			for( int j = 0; j < 4; j++ )
				if( indices[ j ] == ele )
					return j;
			return -1;
		}

		inline int find_opp( int v0_id, int v1_id, int v2_id ) const
		{
			for( int j = 0; j < 4; j++ )
				if( indices[ j ] != v0_id && indices[ j ] != v1_id && indices[ j ] != v2_id )
					return j;
			return -1;
		}

		inline void print( const Logger& log, int id ) const
		{
			log.logDebug( "%i: [ %i, %i, %i, %i ]", id, indices[ 0 ], indices[ 1 ], indices[ 2 ], indices[ 3 ] );
		}

		std::array<int8_t, 4> is_surface_fs = { { NOT_SURFACE, NOT_SURFACE, NOT_SURFACE, NOT_SURFACE } };
		std::array<char, 4> is_bbox_fs = { { NOT_BBOX, NOT_BBOX, NOT_BBOX, NOT_BBOX } };
		std::array<int, 4> opp_t_ids = { { OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN, OPP_T_ID_UNKNOWN } };
		std::array<char, 4> surface_tags = { { 0, 0, 0, 0 } };

		Scalar quality = 0;
		Scalar scalar = 0;
		bool is_removed = false;
		bool is_outside = false;
	};
}