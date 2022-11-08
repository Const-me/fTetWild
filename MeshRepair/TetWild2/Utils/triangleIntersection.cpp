#include <stdafx.h>
#include "triangleIntersection.h"
#include "Geogram2.h"

extern "C++" int tri_tri_intersection_test_3d(
  double p1[ 3 ], double q1[ 3 ], double r1[ 3 ], double p2[ 3 ], double q2[ 3 ], double r2[ 3 ], int* coplanar, double source[ 3 ], double target[ 3 ] );

namespace
{
	int triangleIntersectionTestV1( const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2, bool* coplanar,
	  double* source, double* target )
	{
		int cp;
		int res = tri_tri_intersection_test_3d( const_cast<double*>( p1 ), const_cast<double*>( q1 ), const_cast<double*>( r1 ), const_cast<double*>( p2 ),
		  const_cast<double*>( q2 ), const_cast<double*>( r2 ), &cp, source, target );
		*coplanar = ( 0 != cp );
		return res;
	}

#define SUB( dest, v1, v2 )        \
	dest[ 0 ] = v1[ 0 ] - v2[ 0 ]; \
	dest[ 1 ] = v1[ 1 ] - v2[ 1 ]; \
	dest[ 2 ] = v1[ 2 ] - v2[ 2 ]

#define CROSS( dest, v1, v2 )                          \
	dest[ 0 ] = v1[ 1 ] * v2[ 2 ] - v1[ 2 ] * v2[ 1 ]; \
	dest[ 1 ] = v1[ 2 ] * v2[ 0 ] - v1[ 0 ] * v2[ 2 ]; \
	dest[ 2 ] = v1[ 0 ] * v2[ 1 ] - v1[ 1 ] * v2[ 0 ]

#define DOT( v1, v2 ) ( v1[ 0 ] * v2[ 0 ] + v1[ 1 ] * v2[ 1 ] + v1[ 2 ] * v2[ 2 ] )

#define SCALAR( dest, alpha, v ) \
	dest[ 0 ] = alpha * v[ 0 ];  \
	dest[ 1 ] = alpha * v[ 1 ];  \
	dest[ 2 ] = alpha * v[ 2 ]

	inline int sub_sub_cross_sub_dot( const double* pa, const double* pb, const double* pc, const double* pd )
	{
		auto result = -GEO2::orient_3d( pa, pb, pc, pd );
		if( result > 0 )
			return 1;
		else if( result < 0 )
			return -1;
		return 0;
	}

	template<int dp1, int dq1, int dr1>
	__forceinline uint8_t fsc()
	{
		if( dp1 > 0 )
		{
			if( dq1 > 0 )
				// TRI_TRI_INTER_3D( r1, p1, q1, p2, r2, q2, dp2, dr2, dq2 )
				return 0;
			else if( dr1 > 0 )
				// TRI_TRI_INTER_3D( q1, r1, p1, p2, r2, q2, dp2, dr2, dq2 )
				return 1;
			else
				// TRI_TRI_INTER_3D( p1, q1, r1, p2, q2, r2, dp2, dq2, dr2 )
				return 2;
		}
		else if( dp1 < 0 )
		{
			if( dq1 < 0 )
				// TRI_TRI_INTER_3D( r1, p1, q1, p2, q2, r2, dp2, dq2, dr2 )
				return 3;
			else if( dr1 < 0 )
				// TRI_TRI_INTER_3D( q1, r1, p1, p2, q2, r2, dp2, dq2, dr2 )
				return 4;
			else
				// TRI_TRI_INTER_3D( p1, q1, r1, p2, r2, q2, dp2, dr2, dq2 )
				return 5;
		}
		else
		{
			if( dq1 < 0 )
			{
				if( dr1 >= 0 )
					// TRI_TRI_INTER_3D( q1, r1, p1, p2, r2, q2, dp2, dr2, dq2 )
					return 1;
				else
					// TRI_TRI_INTER_3D( p1, q1, r1, p2, q2, r2, dp2, dq2, dr2 )
					return 2;
			}
			else if( dq1 > 0 )
			{
				if( dr1 > 0 )
					// TRI_TRI_INTER_3D( p1, q1, r1, p2, r2, q2, dp2, dr2, dq2 )
					return 5;
				else
					// TRI_TRI_INTER_3D( q1, r1, p1, p2, q2, r2, dp2, dq2, dr2 )
					return 4;
			}
			else
			{
				if( dr1 > 0 )
					// TRI_TRI_INTER_3D( r1, p1, q1, p2, q2, r2, dp2, dq2, dr2 )
					return 3;
				else if( dr1 < 0 )
					// TRI_TRI_INTER_3D( r1, p1, q1, p2, r2, q2, dp2, dr2, dq2 )
					return 0;
				else
				{
					// triangles are co-planar
					return 0xFF;
				}
			}
		}
	}

	// clang-format off
	static const std::array<uint8_t,3*3*3> g_lookupOuter = 
	{
	  fsc<-1, -1, -1>(), fsc<-1, -1, 0>(), fsc<-1, -1, +1>(),
	  fsc<-1,  0, -1>(), fsc<-1,  0, 0>(), fsc<-1,  0, +1>(),
	  fsc<-1, +1, -1>(), fsc<-1, +1, 0>(), fsc<-1, +1, +1>(),

	  fsc< 0, -1, -1>(), fsc< 0, -1, 0>(), fsc< 0, -1, +1>(),
	  fsc< 0,  0, -1>(), fsc< 0,  0, 0>(), fsc< 0,  0, +1>(),
	  fsc< 0, +1, -1>(), fsc< 0, +1, 0>(), fsc< 0, +1, +1>(),

	  fsc<+1, -1, -1>(), fsc<+1, -1, 0>(), fsc<+1, -1, +1>(),
	  fsc<+1,  0, -1>(), fsc<+1,  0, 0>(), fsc<+1,  0, +1>(),
	  fsc<+1, +1, -1>(), fsc<+1, +1, 0>(), fsc<+1, +1, +1>(),
	};

	inline uint8_t lookupShuffleIndex( int p, int q, int r )
	{
		assert( p >= -1 && p <= +1 );
		assert( q >= -1 && q <= +1 );
		assert( r >= -1 && r <= +1 );

		const size_t idx = (size_t)( p + 1 ) * 9 + (size_t)( q + 1 ) * 3 + (size_t)( r + 1 );
		return g_lookupOuter[ idx ];
	}

	__forceinline std::array<uint8_t, 6> shuffle( uint8_t e0, uint8_t e1, uint8_t e2, uint8_t e3, uint8_t e4, uint8_t e5 )
	{
		return std::array<uint8_t, 6>{ e0, e1, e2, e3, e4, e5 };
	};

	static const std::array<std::array<uint8_t, 6>, 6> switch1shuffles =
	{
		// TRI_TRI_INTER_3D( r1, p1, q1, p2, r2, q2, dp2, dr2, dq2 )
		shuffle( 2, 0, 1,  0, 2, 1  ),
		// TRI_TRI_INTER_3D( q1, r1, p1, p2, r2, q2, dp2, dr2, dq2 )
		shuffle( 1, 2, 0,  0, 2, 1  ),
		// TRI_TRI_INTER_3D( p1, q1, r1, p2, q2, r2, dp2, dq2, dr2 )
		shuffle( 0, 1, 2,  0, 1, 2  ),

		// TRI_TRI_INTER_3D( r1, p1, q1, p2, q2, r2, dp2, dq2, dr2 )
		shuffle( 2, 0, 1,  0, 1, 2  ),
		// TRI_TRI_INTER_3D( q1, r1, p1, p2, q2, r2, dp2, dq2, dr2 )
		shuffle( 1, 2, 0,  0, 1, 2  ),
		// TRI_TRI_INTER_3D( p1, q1, r1, p2, r2, q2, dp2, dr2, dq2 )
		shuffle( 0, 1, 2,  0, 2, 1  ),
/*
		// TRI_TRI_INTER_3D( q1, r1, p1, p2, r2, q2, dp2, dr2, dq2 ) - duplicate of #1
		shuffle( 1, 2, 0,  0, 2, 1  ),
		// TRI_TRI_INTER_3D( p1, q1, r1, p2, q2, r2, dp2, dq2, dr2 ) - duplicate of #2
		shuffle( 0, 1, 2,  0, 1, 2  ),
		// TRI_TRI_INTER_3D( p1, q1, r1, p2, r2, q2, dp2, dr2, dq2 ) - duplicate of #5
		shuffle( 0, 1, 2,  0, 2, 1  ),
		// TRI_TRI_INTER_3D( q1, r1, p1, p2, q2, r2, dp2, dq2, dr2 ) - duplicate of #4
		shuffle( 1, 2, 0,  0, 1, 2  ),
		// TRI_TRI_INTER_3D( r1, p1, q1, p2, q2, r2, dp2, dq2, dr2 ) - duplicate of #3
		shuffle( 2, 0, 1,  0, 1, 2  ),
		// TRI_TRI_INTER_3D( r1, p1, q1, p2, r2, q2, dp2, dr2, dq2 ) - duplicate of #0
		shuffle( 2, 0, 1,  0, 2, 1  ),
*/
	};

	static const std::array<std::array<uint8_t, 6>, 6> switch2shuffles
	{
		// CONSTRUCT_INTERSECTION( p1, r1, q1, r2, p2, q2 )
		shuffle( 0, 2, 1,  2, 0, 1  ),
		// CONSTRUCT_INTERSECTION( p1, r1, q1, q2, r2, p2 )
		shuffle( 0, 2, 1,  1, 2, 0  ),
		// CONSTRUCT_INTERSECTION( p1, q1, r1, p2, q2, r2 )
		shuffle( 0, 1, 2,  0, 1, 2  ),
		// CONSTRUCT_INTERSECTION( p1, q1, r1, r2, p2, q2 )
		shuffle( 0, 1, 2,  2, 0, 1  ),
		// CONSTRUCT_INTERSECTION( p1, q1, r1, q2, r2, p2 )
		shuffle( 0, 1, 2,  1, 2, 0  ),
		// CONSTRUCT_INTERSECTION( p1, r1, q1, p2, q2, r2 )
		shuffle( 0, 2, 1,  0, 1, 2  ),
	};

	// clang-format on

	// CONSTRUCT_INTERSECTION macro, reworked into a function
	int constructIntersection( const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2,
	  const std::array<double, 3>& N1, const std::array<double, 3>& N2, double* source, double* target )
	{
		std::array<double, 3> v1;
		std::array<double, 3> v2;
		double alpha;

		if( sub_sub_cross_sub_dot( q1, r2, p1, p2 ) > 0 )
		{
			if( sub_sub_cross_sub_dot( r1, r2, p1, p2 ) <= 0 )
			{
				if( sub_sub_cross_sub_dot( r1, q2, p1, p2 ) > 0 )
				{
					SUB( v1, p1, p2 );
					SUB( v2, p1, r1 );
					alpha = DOT( v1, N2 ) / DOT( v2, N2 );
					SCALAR( v1, alpha, v2 );
					SUB( source, p1, v1 );
					SUB( v1, p2, p1 );
					SUB( v2, p2, r2 );
					alpha = DOT( v1, N1 ) / DOT( v2, N1 );
					SCALAR( v1, alpha, v2 );
					SUB( target, p2, v1 );
					return 1;
				}
				else
				{
					SUB( v1, p2, p1 );
					SUB( v2, p2, q2 );
					alpha = DOT( v1, N1 ) / DOT( v2, N1 );
					SCALAR( v1, alpha, v2 );
					SUB( source, p2, v1 );
					SUB( v1, p2, p1 );
					SUB( v2, p2, r2 );
					alpha = DOT( v1, N1 ) / DOT( v2, N1 );
					SCALAR( v1, alpha, v2 );
					SUB( target, p2, v1 );
					return 1;
				}
			}
			else
				return 0;
		}
		else
		{
			if( sub_sub_cross_sub_dot( q1, q2, p1, p2 ) < 0 )
				return 0;

			if( sub_sub_cross_sub_dot( r1, q2, p1, p2 ) >= 0 )
			{
				SUB( v1, p1, p2 );
				SUB( v2, p1, r1 );
				alpha = DOT( v1, N2 ) / DOT( v2, N2 );
				SCALAR( v1, alpha, v2 );
				SUB( source, p1, v1 );
				SUB( v1, p1, p2 );
				SUB( v2, p1, q1 );
				alpha = DOT( v1, N2 ) / DOT( v2, N2 );
				SCALAR( v1, alpha, v2 );
				SUB( target, p1, v1 );
				return 1;
			}
			else
			{
				SUB( v1, p2, p1 );
				SUB( v2, p2, q2 );
				alpha = DOT( v1, N1 ) / DOT( v2, N1 );
				SCALAR( v1, alpha, v2 );
				SUB( source, p2, v1 );
				SUB( v1, p1, p2 );
				SUB( v2, p1, q1 );
				alpha = DOT( v1, N2 ) / DOT( v2, N2 );
				SCALAR( v1, alpha, v2 );
				SUB( target, p1, v1 );
				return 1;
			}
		}
	}

	// TRI_TRI_INTER_3D macro, reworked into a function
	int triInter3D( const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2, int dp2, int dq2, int dr2,
	  const std::array<double, 3>& N1, const std::array<double, 3>& N2, bool* coplanar, double* source, double* target )
	{
		const uint8_t idx = lookupShuffleIndex( dp2, dq2, dr2 );
		if( idx == 0xFF )
		{
			// triangles are co-planar
			*coplanar = true;
			return -1;
		}

		const std::array<uint8_t, 6>& shuff = switch2shuffles[ idx ];
		const std::array<const double*, 3> pqr1 = { p1, q1, r1 };
		const std::array<const double*, 3> pqr2 = { p2, q2, r2 };

		return constructIntersection(
		  pqr1[ shuff[ 0 ] ], pqr1[ shuff[ 1 ] ], pqr1[ shuff[ 2 ] ], pqr2[ shuff[ 3 ] ], pqr2[ shuff[ 4 ] ], pqr2[ shuff[ 5 ] ], N1, N2, source, target );
	}

	int triangleIntersectionTestV2( const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2, bool* coplanar,
	  double* source, double* target )
	{
		std::array<int, 3> dpqr1, dpqr2;
		std::array<double, 3> v1, v2;
		std::array<double, 3> N1, N2;

		SUB( v1, q1, p1 );
		SUB( v2, r1, p1 );
		CROSS( N1, v1, v2 );

		SUB( v1, p2, r2 );
		SUB( v2, q2, r2 );
		CROSS( N2, v1, v2 );

		dpqr1[ 0 ] = sub_sub_cross_sub_dot( p2, q2, r2, p1 );
		dpqr1[ 1 ] = sub_sub_cross_sub_dot( p2, q2, r2, q1 );
		dpqr1[ 2 ] = sub_sub_cross_sub_dot( p2, q2, r2, r1 );

		if( ( ( dpqr1[ 0 ] * dpqr1[ 1 ] ) > 0 ) && ( ( dpqr1[ 0 ] * dpqr1[ 2 ] ) > 0 ) )
			return 666;

		// Compute distance signs  of p2, q2 and r2
		// to the plane of triangle(p1,q1,r1)

		dpqr2[ 0 ] = sub_sub_cross_sub_dot( p1, q1, r1, p2 );
		dpqr2[ 1 ] = sub_sub_cross_sub_dot( p1, q1, r1, q2 );
		dpqr2[ 2 ] = sub_sub_cross_sub_dot( p1, q1, r1, r2 );
		if( ( ( dpqr2[ 0 ] * dpqr2[ 1 ] ) > 0 ) && ( ( dpqr2[ 0 ] * dpqr2[ 2 ] ) > 0 ) )
			return 666;

		// The first bunch of IFs
		const uint8_t index2 = lookupShuffleIndex( dpqr1[ 0 ], dpqr1[ 1 ], dpqr1[ 2 ] );
		if( index2 == 0xFF )
		{
			// triangles are co-planar
			*coplanar = true;
			return -1;
		}
		const std::array<uint8_t, 6>& shuff = switch1shuffles[ index2 ];
		const std::array<const double*, 3> pqr1 = { p1, q1, r1 };
		const std::array<const double*, 3> pqr2 = { p2, q2, r2 };
		return triInter3D( pqr1[ shuff[ 0 ] ], pqr1[ shuff[ 1 ] ], pqr1[ shuff[ 2 ] ], pqr2[ shuff[ 3 ] ], pqr2[ shuff[ 4 ] ], pqr2[ shuff[ 5 ] ],
		  dpqr2[ shuff[ 3 ] ], dpqr2[ shuff[ 4 ] ], dpqr2[ shuff[ 5 ] ], N1, N2, coplanar, source, target );
	}

	inline bool equals( const std::array<double, 3>& a, const double* p )
	{
		return 0 == memcmp( a.data(), p, 8 * 3 );
	}
}  // namespace

int triangleIntersectionTest(
  const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2, bool* coplanar, double* source, double* target )
{
#if 1
	return triangleIntersectionTestV2( p1, q1, r1, p2, q2, r2, coplanar, source, target );
#else
	std::array<double, 3> s2 { source[ 0 ], source[ 1 ], source[ 2 ] };
	std::array<double, 3> t2 { target[ 0 ], target[ 1 ], target[ 2 ] };
	bool cp2 = *coplanar;

	const int res2 = triangleIntersectionTestV2( p1, q1, r1, p2, q2, r2, &cp2, s2.data(), t2.data() );
	const int res1 = triangleIntersectionTestV1( p1, q1, r1, p2, q2, r2, coplanar, source, target );

	int diff = 0;
	if( res1 != res2 )
		diff |= 1;
	if( *coplanar != cp2 )
		diff |= 2;
	if( !equals( s2, source ) )
		diff |= 4;
	if( !equals( t2, target ) )
		diff |= 8;
	if( 0 == diff )
		return res2;

	__debugbreak();
	triangleIntersectionTestV2( p1, q1, r1, p2, q2, r2, &cp2, s2.data(), t2.data() );
	return res2;
#endif
}