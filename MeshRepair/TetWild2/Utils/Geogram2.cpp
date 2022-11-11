#include "stdafx.h"
#include "Geogram2.h"
#include <geogram/basic/geometry_nd.h>
#ifdef __AVX__
#include "AvxMath.h"
#else
#error This code requires at least AVX1
#endif

namespace GEO2
{
	static double pointTriangleSquaredDistanceOrig( const vec3& point, const vec3& V0, const vec3& V1, const vec3& V2, vec3* closest_point )
	{
		vec3 cp;
		double lambda0, lambda1, lambda2;
		double res = GEO::Geom::point_triangle_squared_distance( point, V0, V1, V2, cp, lambda0, lambda1, lambda2 );
		if( nullptr != closest_point )
			*closest_point = cp;
		return res;
	}

	namespace
	{
		__forceinline double dot( __m256d a, __m256d b )
		{
			return AvxMath::vector3DotScalar( a, b );
		}
		__forceinline double length2( __m256d v )
		{
			return AvxMath::vector3DotScalar( v, v );
		}
		__forceinline double distance2( __m256d a, __m256d b )
		{
			__m256d v = _mm256_sub_pd( b, a );
			return length2( v );
		}
	}  // namespace

	inline double pointSegmentSquaredDistance( const __m256d point, const __m256d v0, const __m256d v1, __m256d& closest )
	{
		double l2 = distance2( v0, v1 );
		double t = dot( _mm256_sub_pd( point, v0 ), _mm256_sub_pd( v1, v0 ) );
		if( t <= 0.0 || l2 == 0.0 )
		{
			closest = v0;
			return distance2( point, v0 );
		}
		else if( t > l2 )
		{
			closest = v1;
			return distance2( point, v1 );
		}

		__m256d cp = AvxMath::lerpFast( v0, v1, t / l2 );
		closest = cp;
		return distance2( point, cp );
	}

	static double __declspec( noinline )
	  pointDegenerateTriangleSquaredDistance( const __m256d point, const __m256d v0, const __m256d v1, const __m256d v2, vec3* closestPoint )
	{
		__m256d closest;
		double result = pointSegmentSquaredDistance( point, v0, v1, closest );

		__m256d cp2;
		double r2 = pointSegmentSquaredDistance( point, v0, v2, cp2 );
		if( r2 < result )
		{
			result = r2;
			closest = cp2;
		}

		r2 = pointSegmentSquaredDistance( point, v1, v2, cp2 );
		if( r2 < result )
		{
			result = r2;
			closest = cp2;
		}

		if( nullptr != closestPoint )
			AvxMath::storeDouble3( &closestPoint->x, closest );

		return result;
	}

	double pointTriangleSquaredDistanceAvx( const vec3& Point, const vec3& V0, const vec3& V1, const vec3& V2, vec3* closest_point )
	{
		using namespace AvxMath;
		const __m256d point = loadDouble3( &Point.x );
		const __m256d v0 = loadDouble3( &V0.x );
		const __m256d v1 = loadDouble3( &V1.x );
		const __m256d v2 = loadDouble3( &V2.x );

		const __m256d diff = _mm256_sub_pd( v0, point );
		const __m256d edge0 = _mm256_sub_pd( v1, v0 );
		const __m256d edge1 = _mm256_sub_pd( v2, v0 );

		const double a00 = length2( edge0 );
		const double a01 = dot( edge0, edge1 );
		const double a11 = length2( edge1 );
		const double b0 = dot( diff, edge0 );
		const double b1 = dot( diff, edge1 );
		const double c = length2( diff );
		const double det = ::fabs( a00 * a11 - a01 * a01 );

		// If the triangle is degenerate
		if( det < 1e-30 )
			return pointDegenerateTriangleSquaredDistance( point, v0, v1, v2, closest_point );

		double s = a01 * b1 - a11 * b0;
		double t = a01 * b0 - a00 * b1;
		double sqrDistance;
		if( s + t <= det )
		{
			if( s < 0.0 )
			{
				if( t < 0.0 )
				{  // region 4
					if( b0 < 0.0 )
					{
						t = 0.0;
						if( -b0 >= a00 )
						{
							s = 1.0;
							sqrDistance = a00 + 2.0 * b0 + c;
						}
						else
						{
							s = -b0 / a00;
							sqrDistance = b0 * s + c;
						}
					}
					else
					{
						s = 0.0;
						if( b1 >= 0.0 )
						{
							t = 0.0;
							sqrDistance = c;
						}
						else if( -b1 >= a11 )
						{
							t = 1.0;
							sqrDistance = a11 + 2.0 * b1 + c;
						}
						else
						{
							t = -b1 / a11;
							sqrDistance = b1 * t + c;
						}
					}
				}
				else
				{  // region 3
					s = 0.0;
					if( b1 >= 0.0 )
					{
						t = 0.0;
						sqrDistance = c;
					}
					else if( -b1 >= a11 )
					{
						t = 1.0;
						sqrDistance = a11 + 2.0 * b1 + c;
					}
					else
					{
						t = -b1 / a11;
						sqrDistance = b1 * t + c;
					}
				}
			}
			else if( t < 0.0 )
			{  // region 5
				t = 0.0;
				if( b0 >= 0.0 )
				{
					s = 0.0;
					sqrDistance = c;
				}
				else if( -b0 >= a00 )
				{
					s = 1.0;
					sqrDistance = a00 + 2.0 * b0 + c;
				}
				else
				{
					s = -b0 / a00;
					sqrDistance = b0 * s + c;
				}
			}
			else
			{  // region 0
				// minimum at interior point
				double invDet = double( 1.0 ) / det;
				s *= invDet;
				t *= invDet;
				sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 ) + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c;
			}
		}
		else
		{
			double tmp0, tmp1, numer, denom;

			if( s < 0.0 )
			{  // region 2
				tmp0 = a01 + b0;
				tmp1 = a11 + b1;
				if( tmp1 > tmp0 )
				{
					numer = tmp1 - tmp0;
					denom = a00 - 2.0 * a01 + a11;
					if( numer >= denom )
					{
						s = 1.0;
						t = 0.0;
						sqrDistance = a00 + 2.0 * b0 + c;
					}
					else
					{
						s = numer / denom;
						t = 1.0 - s;
						sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 ) + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c;
					}
				}
				else
				{
					s = 0.0;
					if( tmp1 <= 0.0 )
					{
						t = 1.0;
						sqrDistance = a11 + 2.0 * b1 + c;
					}
					else if( b1 >= 0.0 )
					{
						t = 0.0;
						sqrDistance = c;
					}
					else
					{
						t = -b1 / a11;
						sqrDistance = b1 * t + c;
					}
				}
			}
			else if( t < 0.0 )
			{  // region 6
				tmp0 = a01 + b1;
				tmp1 = a00 + b0;
				if( tmp1 > tmp0 )
				{
					numer = tmp1 - tmp0;
					denom = a00 - 2.0 * a01 + a11;
					if( numer >= denom )
					{
						t = 1.0;
						s = 0.0;
						sqrDistance = a11 + 2.0 * b1 + c;
					}
					else
					{
						t = numer / denom;
						s = 1.0 - t;
						sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 ) + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c;
					}
				}
				else
				{
					t = 0.0;
					if( tmp1 <= 0.0 )
					{
						s = 1.0;
						sqrDistance = a00 + 2.0 * b0 + c;
					}
					else if( b0 >= 0.0 )
					{
						s = 0.0;
						sqrDistance = c;
					}
					else
					{
						s = -b0 / a00;
						sqrDistance = b0 * s + c;
					}
				}
			}
			else
			{  // region 1
				numer = a11 + b1 - a01 - b0;
				if( numer <= 0.0 )
				{
					s = 0.0;
					t = 1.0;
					sqrDistance = a11 + 2.0 * b1 + c;
				}
				else
				{
					denom = a00 - 2.0 * a01 + a11;
					if( numer >= denom )
					{
						s = 1.0;
						t = 0.0;
						sqrDistance = a00 + 2.0 * b0 + c;
					}
					else
					{
						s = numer / denom;
						t = 1.0 - s;
						sqrDistance = s * ( a00 * s + a01 * t + 2.0 * b0 ) + t * ( a01 * s + a11 * t + 2.0 * b1 ) + c;
					}
				}
			}
		}
		// Account for numerical round-off error.
		if( sqrDistance < 0.0 )
			sqrDistance = 0.0;

		if( nullptr != closest_point )
		{
			__m256d tmp1 = _mm256_mul_pd( edge0, _mm256_set1_pd( s ) );
			__m256d tmp2 = _mm256_mul_pd( edge1, _mm256_set1_pd( t ) );
			__m256d cp = _mm256_add_pd( _mm256_add_pd( v0, tmp1 ), tmp2 );
			storeDouble3( &closest_point->x, cp );
		}

		return sqrDistance;
	}

	double point_triangle_squared_distance( const vec3& point, const vec3& V0, const vec3& V1, const vec3& V2, vec3* closest_point )
	{
#if 1
		return pointTriangleSquaredDistanceAvx( point, V0, V1, V2, closest_point );
#else
		double orig = pointTriangleSquaredDistanceOrig( point, V0, V1, V2, closest_point );
		double my = pointTriangleSquaredDistanceAvx( point, V0, V1, V2, closest_point );
		if( my != orig )
			__debugbreak();
		return orig;
#endif
	}

	static void pointTriDist()
	{
		vec3 pt { 8.2750787734985352, 34.215789794921875, 43.730640411376953 };
		vec3 v0 { 8.2715597152709961, 34.215789794921875, 43.764884948730469 };
		vec3 v1 { 8.2715597152709961, 34.215789794921875, 43.764884948730469 };
		vec3 v2 { 8.2750787734985352, 34.215789794921875, 43.730640411376953 };
		point_triangle_squared_distance( pt, v0, v1, v2 );
	}

	void dbgRunSomeTests()
	{
		pointTriDist();
	}
}  // namespace GEO2