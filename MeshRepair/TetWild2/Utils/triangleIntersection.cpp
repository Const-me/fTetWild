// Fast and Robust Triangle-Triangle Overlap Test Using Orientation Predicates
// Philippe Guigue, Olivier Devillers, January 2003
// DOI: 10.1080/10867651.2003.10487580
#include <stdafx.h>
#include "triangleIntersection.h"
#include "Geogram2.h"
#include "SimdVector.h"

namespace
{
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
				// CONSTRUCT_INTERSECTION( p1, r1, q1, r2, p2, q2 )
				return 0;
			else if( dr1 > 0 )
				// TRI_TRI_INTER_3D( q1, r1, p1, p2, r2, q2, dp2, dr2, dq2 )
				// CONSTRUCT_INTERSECTION( p1, r1, q1, q2, r2, p2 )
				return 1;
			else
				// TRI_TRI_INTER_3D( p1, q1, r1, p2, q2, r2, dp2, dq2, dr2 )
				// CONSTRUCT_INTERSECTION( p1, q1, r1, p2, q2, r2 )
				return 2;
		}
		else if( dp1 < 0 )
		{
			if( dq1 < 0 )
				// TRI_TRI_INTER_3D( r1, p1, q1, p2, q2, r2, dp2, dq2, dr2 )
				// CONSTRUCT_INTERSECTION( p1, q1, r1, r2, p2, q2 )
				return 3;
			else if( dr1 < 0 )
				// TRI_TRI_INTER_3D( q1, r1, p1, p2, q2, r2, dp2, dq2, dr2 )
				// CONSTRUCT_INTERSECTION( p1, q1, r1, q2, r2, p2 )
				return 4;
			else
				// TRI_TRI_INTER_3D( p1, q1, r1, p2, r2, q2, dp2, dr2, dq2 )
				// CONSTRUCT_INTERSECTION( p1, r1, q1, p2, q2, r2 )
				return 5;
		}
		else
		{
			if( dq1 < 0 )
			{
				if( dr1 >= 0 )
					// TRI_TRI_INTER_3D( q1, r1, p1, p2, r2, q2, dp2, dr2, dq2 )
					// CONSTRUCT_INTERSECTION( p1, r1, q1, q2, r2, p2 )
					return 1;
				else
					// TRI_TRI_INTER_3D( p1, q1, r1, p2, q2, r2, dp2, dq2, dr2 )
					// CONSTRUCT_INTERSECTION( p1, q1, r1, p2, q2, r2 )
					return 2;
			}
			else if( dq1 > 0 )
			{
				if( dr1 > 0 )
					// TRI_TRI_INTER_3D( p1, q1, r1, p2, r2, q2, dp2, dr2, dq2 )
					// CONSTRUCT_INTERSECTION( p1, r1, q1, p2, q2, r2 )
					return 5;
				else
					// TRI_TRI_INTER_3D( q1, r1, p1, p2, q2, r2, dp2, dq2, dr2 )
					// CONSTRUCT_INTERSECTION( p1, q1, r1, q2, r2, p2 )
					return 4;
			}
			else
			{
				if( dr1 > 0 )
					// TRI_TRI_INTER_3D( r1, p1, q1, p2, q2, r2, dp2, dq2, dr2 )
					// CONSTRUCT_INTERSECTION( p1, q1, r1, r2, p2, q2 )
					return 3;
				else if( dr1 < 0 )
					// TRI_TRI_INTER_3D( r1, p1, q1, p2, r2, q2, dp2, dr2, dq2 )
					// CONSTRUCT_INTERSECTION( p1, r1, q1, r2, p2, q2 )
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
	using Vec = Simd::Vec3;
	using Simd::cross;
	using Simd::dot;
	using Simd::load3;
	using Simd::store3;

	struct Context
	{
		Vec N1;
		Vec N2;
		double* source;
		double* target;
		bool* coplanar;
	};

	// CONSTRUCT_INTERSECTION macro, reworked into a function
	int constructIntersection( Context& c, const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2 )
	{
		if( sub_sub_cross_sub_dot( q1, r2, p1, p2 ) > 0 )
		{
			if( sub_sub_cross_sub_dot( r1, r2, p1, p2 ) <= 0 )
			{
				if( sub_sub_cross_sub_dot( r1, q2, p1, p2 ) > 0 )
				{
					const Vec p1v = load3( p1 );
					const Vec p2v = load3( p2 );
					Vec v1 = p1v - p2v;
					Vec v2 = p1v - load3( r1 );
					double alpha = dot( v1, c.N2 ) / dot( v2, c.N2 );
					v1 = v2 * alpha;
					store3( c.source, p1v - v1 );
					v1 = p2v - p1v;
					v2 = p2v - load3( r2 );
					alpha = dot( v1, c.N1 ) / dot( v2, c.N1 );
					v1 = v2 * alpha;
					store3( c.target, p2v - v1 );
					return 1;
				}
				else
				{
					const Vec p1v = load3( p1 );
					const Vec p2v = load3( p2 );
					Vec v1 = p2v - p1v;
					Vec v2 = p2v - load3( q2 );
					double alpha = dot( v1, c.N1 ) / dot( v2, c.N1 );
					v1 = v2 * alpha;
					store3( c.source, p2v - v1 );
					v1 = p2v - p1v;
					v2 = p2v - load3( r2 );
					alpha = dot( v1, c.N1 ) / dot( v2, c.N1 );
					v1 = v2 * alpha;
					store3( c.target, p2v - v1 );
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
				const Vec p1v = load3( p1 );
				const Vec p2v = load3( p2 );
				Vec v1 = p1v - p2v;
				Vec v2 = p1v - load3( r1 );
				double alpha = dot( v1, c.N2 ) / dot( v2, c.N2 );
				v1 = v2 * alpha;
				store3( c.source, p1v - v1 );
				v1 = p1v - p2v;
				v2 = p1v - load3( q1 );
				alpha = dot( v1, c.N2 ) / dot( v2, c.N2 );
				v1 = v2 * alpha;
				store3( c.target, p1v - v1 );
				return 1;
			}
			else
			{
				const Vec p1v = load3( p1 );
				const Vec p2v = load3( p2 );

				Vec v1 = p2v - p1v;
				Vec v2 = p2v - load3( q2 );
				double alpha = dot( v1, c.N1 ) / dot( v2, c.N1 );
				v1 = v2 * alpha;
				store3( c.source, p2v - v1 );
				v1 = p1v - p2v;
				v2 = p1v - load3( q1 );
				alpha = dot( v1, c.N2 ) / dot( v2, c.N2 );
				v1 = v2 * alpha;
				store3( c.target, p1v - v1 );
				return 1;
			}
		}
	}

	// TRI_TRI_INTER_3D macro, reworked into a function
	int triInter3D(
	  Context& c, const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2, int dp2, int dq2, int dr2 )
	{
		const uint8_t idx = lookupShuffleIndex( dp2, dq2, dr2 );
		if( idx == 0xFF )
		{
			// triangles are co-planar
			*c.coplanar = true;
			return -1;
		}

		const std::array<uint8_t, 6>& shuff = switch2shuffles[ idx ];
		const std::array<const double*, 3> pqr1 = { p1, q1, r1 };
		const std::array<const double*, 3> pqr2 = { p2, q2, r2 };

		return constructIntersection(
		  c, pqr1[ shuff[ 0 ] ], pqr1[ shuff[ 1 ] ], pqr1[ shuff[ 2 ] ], pqr2[ shuff[ 3 ] ], pqr2[ shuff[ 4 ] ], pqr2[ shuff[ 5 ] ] );
	}
}

int triangleIntersectionTest(
  const double* p1, const double* q1, const double* r1, const double* p2, const double* q2, const double* r2, bool* coplanar,
	double* source, double* target )
{
	std::array<int, 3> dpqr1, dpqr2;

	dpqr1[ 0 ] = sub_sub_cross_sub_dot( p2, q2, r2, p1 );
	dpqr1[ 1 ] = sub_sub_cross_sub_dot( p2, q2, r2, q1 );
	dpqr1[ 2 ] = sub_sub_cross_sub_dot( p2, q2, r2, r1 );

	if( ( ( dpqr1[ 0 ] * dpqr1[ 1 ] ) > 0 ) && ( ( dpqr1[ 0 ] * dpqr1[ 2 ] ) > 0 ) )
		return 666;

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

	Context c;
	c.coplanar = coplanar;
	c.source = source;
	c.target = target;

	Vec tmp = load3( p1 );
	Vec v1 = load3( q1 ) - tmp;
	Vec v2 = load3( r1 ) - tmp;
	c.N1 = cross( v1, v2 );

	tmp = load3( r2 );
	v1 = load3( p2 ) - tmp;
	v2 = load3( q2 ) - tmp;
	c.N2 = cross( v1, v2 );

	return triInter3D( c, pqr1[ shuff[ 0 ] ], pqr1[ shuff[ 1 ] ], pqr1[ shuff[ 2 ] ], pqr2[ shuff[ 3 ] ], pqr2[ shuff[ 4 ] ], pqr2[ shuff[ 5 ] ],
		dpqr2[ shuff[ 3 ] ], dpqr2[ shuff[ 4 ] ], dpqr2[ shuff[ 5 ] ] );
}