#include "stdafx.h"
#include "Predicates.h"
#include <geogram/numerics/predicates.h>

namespace floatTetWild
{
	eOrientation Predicates::orient_3d( const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4 )
	{
		const int result = -GEO::PCK::orient_3d( p1.data(), p2.data(), p3.data(), p4.data() );

		if( result > 0 )
			return eOrientation::Positive;
		else if( result < 0 )
			return eOrientation::Negative;
		else
			return eOrientation::Zero;
	}

	eOrientation Predicates::orient_3d_tolerance( const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p )
	{
		const int result = -GEO::PCK::orient_3d( p1.data(), p2.data(), p3.data(), p.data() );
		if( result == 0 )
			return eOrientation::Zero;

		Vector3 n = ( ( p2 - p3 ).cross( p1 - p3 ) ).normalized();
		Scalar d = std::abs( n.dot( p - p1 ) );
		if( d <= SCALAR_ZERO )
			return eOrientation::Zero;

		if( result > 0 )
			return eOrientation::Positive;
		else
			return eOrientation::Negative;
	}

	Scalar Predicates::orient_3d_volume( const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& p4 )
	{
		const int ori = -GEO::PCK::orient_3d( p1.data(), p2.data(), p3.data(), p4.data() );
		if( ori <= 0 )
			return ori;
		else
			return ( p1 - p4 ).dot( ( p2 - p4 ).cross( p3 - p4 ) ) / 6;
	}

	eOrientation Predicates::orient_2d( const Vector2& p1, const Vector2& p2, const Vector2& p3 )
	{
		const int result = -GEO::PCK::orient_2d( p1.data(), p2.data(), p3.data() );
		if( result > 0 )
			return eOrientation::Positive;
		else if( result < 0 )
			return eOrientation::Negative;
		else
			return eOrientation::Zero;
	}
}  // namespace floatTetWild