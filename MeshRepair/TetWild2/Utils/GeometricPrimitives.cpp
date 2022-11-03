#include "stdafx.h"
#include "GeometricPrimitives.h"

namespace GEO2
{
	double inner_point_box_squared_distance( const vec3& p, const Box& B )
	{
		assert( B.contains( p ) );
		/*
		double result = geo_sqr( p[ 0 ] - B.xyz_min[ 0 ] );
		result = std::min( result, geo_sqr( p[ 0 ] - B.xyz_max[ 0 ] ) );
		for( coord_index_t c = 1; c < 3; ++c )
		{
			result = std::min( result, geo_sqr( p[ c ] - B.xyz_min[ c ] ) );
			result = std::min( result, geo_sqr( p[ c ] - B.xyz_max[ c ] ) );
		}
		return result;
		*/
		__debugbreak();
		return 0;
	}

	double tetra_signed_volume( const vec3& p1, const vec3& p2, const vec3& p3, const vec3& p4 )
	{
		return dot( p2 - p1, cross( p3 - p1, p4 - p1 ) ) / 6.0;
	}

	bool Box::contains( const vec3& pos ) const
	{
		const __m128d p1 = _mm_loadu_pd( &pos.x );
		const __m128d p2 = _mm_load_sd( &pos.z );

		__m128d v = _mm_load_pd( &xyz_min[ 0 ] );
		__m128d oob = _mm_cmplt_pd( p1, v );  // p.xy < min.xy
		v = _mm_load_sd( &xyz_min[ 2 ] );
		oob = _mm_or_pd( oob, _mm_cmplt_pd( p2, v ) );	// |= p.z < min.z

		v = _mm_loadu_pd( &xyz_max[ 0 ] );
		oob = _mm_or_pd( oob, _mm_cmpgt_pd( p1, v ) );	// |= p.xy > max.xy
		v = _mm_load_sd( &xyz_max[ 2 ] );
		oob = _mm_or_pd( oob, _mm_cmpgt_pd( p2, v ) );	// |= p.z > max.z

		return (bool)_mm_testz_pd( oob, oob );
	}

	bool bboxes_overlap( const Box& a, const Box& b )
	{
		// a.max.xy < b.min.xy
		__m128d va = _mm_loadu_pd( &a.xyz_max[ 0 ] );
		__m128d vb = _mm_loadu_pd( &b.xyz_min[ 0 ] );
		__m128d oob = _mm_cmplt_pd( va, vb );

		// |= a.max.z < b.min.z
		va = _mm_load_sd( &a.xyz_max[ 2 ] );
		vb = _mm_load_sd( &b.xyz_min[ 2 ] );
		oob = _mm_or_pd( oob, _mm_cmplt_pd( va, vb ) );

		// |= a.min.xy > b.max.xy
		va = _mm_loadu_pd( &a.xyz_min[ 0 ] );
		vb = _mm_loadu_pd( &b.xyz_max[ 0 ] );
		oob = _mm_or_pd( oob, _mm_cmpgt_pd( va, vb ) );

		// |= a.min.z > b.max.z
		va = _mm_load_sd( &a.xyz_min[ 2 ] );
		vb = _mm_load_sd( &b.xyz_max[ 2 ] );
		oob = _mm_or_pd( oob, _mm_cmpgt_pd( va, vb ) );

		return (bool)_mm_testz_pd( oob, oob );
	}

	void bbox_union( Box& rdi, const Box& a, const Box& b )
	{
		// min.xy
		__m128d v = _mm_load_pd( a.xyz_min.data() );
		v = _mm_min_pd( v, _mm_load_pd( b.xyz_min.data() ) );
		_mm_store_pd( rdi.xyz_min.data(), v );
		// min.z
		rdi.xyz_min[ 2 ] = std::min( a.xyz_min[ 2 ], b.xyz_min[ 2 ] );

		// max.x
		rdi.xyz_max[ 0 ] = std::max( a.xyz_max[ 0 ], b.xyz_max[ 0 ] );
		// max.yz - note it's aligned by 16 bytes just like min.xy
		v = _mm_load_pd( &a.xyz_max[ 1 ] );
		v = _mm_max_pd( v, _mm_load_pd( &b.xyz_max[ 1 ] ) );
		_mm_store_pd( &rdi.xyz_max[ 1 ], v );
	}

	static inline double point_segment_squared_distance( const vec3& point, const vec3& V0, const vec3& V1, vec3& closest_point )
	{
		double l2 = distance2( V0, V1 );
		double t = dot( point - V0, V1 - V0 );
		if( t <= 0.0 || l2 == 0.0 )
		{
			closest_point = V0;
			return distance2( point, V0 );
		}
		else if( t > l2 )
		{
			closest_point = V1;
			return distance2( point, V1 );
		}
		double lambda1 = t / l2;
		double lambda0 = 1.0 - lambda1;
		closest_point = lambda0 * V0 + lambda1 * V1;
		return distance2( point, closest_point );
	}

	double point_triangle_squared_distance( const vec3& point, const vec3& V0, const vec3& V1, const vec3& V2, vec3& closest_point )
	{
		vec3 diff = V0 - point;
		vec3 edge0 = V1 - V0;
		vec3 edge1 = V2 - V0;
		double a00 = length2( edge0 );
		double a01 = dot( edge0, edge1 );
		double a11 = length2( edge1 );
		double b0 = dot( diff, edge0 );
		double b1 = dot( diff, edge1 );
		double c = length2( diff );
		double det = ::fabs( a00 * a11 - a01 * a01 );
		double s = a01 * b1 - a11 * b0;
		double t = a01 * b0 - a00 * b1;
		double sqrDistance;

		// If the triangle is degenerate
		if( det < 1e-30 )
		{
			vec3 cur_closest;
			double result;
			double cur_dist = point_segment_squared_distance( point, V0, V1, cur_closest );
			result = cur_dist;
			closest_point = cur_closest;
			cur_dist = point_segment_squared_distance( point, V0, V2, cur_closest );
			if( cur_dist < result )
			{
				result = cur_dist;
				closest_point = cur_closest;
			}
			cur_dist = point_segment_squared_distance( point, V1, V2, cur_closest );
			if( cur_dist < result )
			{
				result = cur_dist;
				closest_point = cur_closest;
			}
			return result;
		}

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

		closest_point = V0 + s * edge0 + t * edge1;
		return sqrDistance;
	}
}