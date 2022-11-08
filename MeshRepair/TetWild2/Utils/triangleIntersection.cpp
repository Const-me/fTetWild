#include <stdafx.h>
#include "triangleIntersection.h"

extern "C++" int tri_tri_intersection_test_3d(
  double p1[ 3 ], double q1[ 3 ], double r1[ 3 ], double p2[ 3 ], double q2[ 3 ], double r2[ 3 ], int* coplanar, double source[ 3 ], double target[ 3 ] );

namespace
{
	bool triangleIntersectionTestV1( const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2, bool* coplanar,
	  double* source, double* target )
	{
		int cp;
		bool res = tri_tri_intersection_test_3d( const_cast<double*>( p1 ), const_cast<double*>( q1 ), const_cast<double*>( r1 ), const_cast<double*>( p2 ),
		  const_cast<double*>( q2 ), const_cast<double*>( r2 ), &cp, source, target );
		*coplanar = ( 0 != cp );
		return ( 0 != res );
	}
}

bool triangleIntersectionTest(
  const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2, bool* coplanar, double* source, double* target )
{
	return triangleIntersectionTestV1( p1, q1, r1, p2, q2, r2, coplanar, source, target );
}