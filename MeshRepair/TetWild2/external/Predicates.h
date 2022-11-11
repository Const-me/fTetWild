#pragma once
#include "../src/Types.hpp"

namespace floatTetWild
{
	enum struct eOrientation : int8_t
	{
		Positive = 1,
		Zero = 0,
		Negative = -1,
		Unknown = 127,
	};

	namespace Predicates
	{
		eOrientation orient_3d( const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4 );
		eOrientation orient_3d_tolerance( const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4 );
		Scalar orient_3d_volume( const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4 );
		eOrientation orient_2d( const Vector2& p1, const Vector2& p2, const Vector2& p3 );
	};	// namespace Predicates
}  // namespace floatTetWild