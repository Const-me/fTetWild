#pragma once
#include "Types.hpp"

namespace floatTetWild
{
	struct MeshVertex
	{
		MeshVertex( const Vector3& p )
			: pos( p )
		{
		}
		MeshVertex()
		{
		}
		Vector3 pos;

		inline Scalar& operator[]( const int index )
		{
			assert( index >= 0 && index < 3 );
			return pos[ index ];
		}

		inline Scalar operator[]( const int index ) const
		{
			assert( index >= 0 && index < 3 );
			return pos[ index ];
		}

		std::vector<int> conn_tets;

		bool is_on_surface = false;
		bool is_on_boundary = false;
		bool is_on_cut = false;
		int on_boundary_e_id = -1;
		bool is_on_bbox = false;
		bool is_outside = false;

		bool is_removed = false;
		bool is_freezed = false;  // todo

		Scalar sizing_scalar = 1;

		Scalar scalar = 0;
	};
}  // namespace floatTetWild