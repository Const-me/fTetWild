#pragma once
#include "Types.hpp"

namespace floatTetWild
{
	// Conditional moves are only there for 16/32/64 bits: https://www.felixcloutier.com/x86/cmovcc
	enum struct eVertexFlags : uint16_t
	{
		Surface = 1,
		Boundary = 2,
		Cut = 4,
		BoundingBox = 8,
		Outside = 0x10,
		Removed = 0x20,
		Freezed = 0x40,
	};

	inline eVertexFlags operator|( eVertexFlags a, eVertexFlags b )
	{
		return (eVertexFlags)( (uint16_t)a | (uint16_t)b );
	}

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
		uint16_t flags = 0;
		inline bool hasFlag( eVertexFlags f ) const
		{
			return 0 != ( flags & (uint16_t)f );
		}
		inline void setFlag( eVertexFlags f, bool val )
		{
			// Hoping to convince the compiler to emit CMOV instruction
			uint16_t s = flags | (uint16_t)f;
			uint16_t us = flags & ( ~( (uint16_t)f ) );
			flags = val ? s : us;
		}
		inline void setFlag( eVertexFlags f )
		{
			flags |= (uint16_t)f;
		}
		inline void clearFlag( eVertexFlags f )
		{
			flags &= ~(uint16_t)f;
		}

		inline bool isRemoved() const
		{
			return hasFlag( eVertexFlags::Removed );
		}
		inline bool isFreezed() const
		{
			return hasFlag( eVertexFlags::Freezed );
		}
		inline bool isSurface() const
		{
			return hasFlag( eVertexFlags::Surface );
		}
		inline bool isCut() const
		{
			return hasFlag( eVertexFlags::Cut );
		}
		inline bool isBoundary() const
		{
			return hasFlag( eVertexFlags::Boundary );
		}
		inline bool isBoundingBox() const
		{
			return hasFlag( eVertexFlags::BoundingBox );
		}

		int on_boundary_e_id = -1;
		Scalar sizing_scalar = 1;
		Scalar scalar = 0;
	};
}  // namespace floatTetWild