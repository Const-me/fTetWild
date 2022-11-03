#pragma once
#include <cmath>

namespace GEO2
{
	struct vec3
	{
		double x, y, z;

		double operator[]( size_t idx ) const
		{
			assert( idx < 3 );
			return ( &x )[ idx ];
		}

		vec3 operator+( const vec3& v ) const
		{
			return vec3 { x + v.x, y + v.y, z + v.z };
		}

		vec3 operator-( const vec3& v ) const
		{
			return vec3 { x - v.x, y - v.y, z - v.z };
		}
		vec3 operator*( double s ) const
		{
			return vec3 { x * s, y * s, z * s };
		}
		vec3 operator/( double s ) const
		{
			return vec3 { x / s, y / s, z / s };
		}

		vec3() = default;
		vec3( double a, double b, double c )
		{
			x = a;
			y = b;
			z = c;
		}
	};

	inline vec3 operator*( double s, const vec3& v )
	{
		return v * s;
	}

	inline double dot( const vec3& a, const vec3& b )
	{
		__m128d av = _mm_loadu_pd( &a.x );
		__m128d bv = _mm_loadu_pd( &b.x );
		__m128d res = _mm_dp_pd( av, bv, 0b00110001 );

		av = _mm_load_sd( &a.z );
		bv = _mm_load_sd( &b.z );
		res = _mm_add_sd( res, _mm_mul_sd( av, bv ) );
		return _mm_cvtsd_f64( res );
	}

	inline double length2( const vec3& vec )
	{
		return dot( vec, vec );
	}
	inline double distance2( const vec3& a, const vec3& b )
	{
		return length2( a - b );
	}
	inline double distance( const vec3& a, const vec3& b )
	{
		return std::sqrt( distance2( a, b ) );
	}
	inline vec3 cross( const vec3& v1, const vec3& v2 )
	{
		return vec3 { v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x };
	}
	inline double length( const vec3& vec )
	{
		return std::sqrt( length2( vec ) );
	}

	inline vec3 normalize( const vec3& v )
	{
		const double s = length( v );
		if( s > 1e-30 )
			return v / s;
		return v * s;
	}

	enum Sign : int
	{
		/** Value is negative */
		NEGATIVE = -1,
		/** Value is zero */
		ZERO = 0,
		/** Value is positive */
		POSITIVE = 1
	};

	inline Sign geo_sgn( double g )
	{
		if( g > 0 )
			return Sign::POSITIVE;
		if( g > 0 )
			return Sign::NEGATIVE;
		return Sign::ZERO;
	}

	struct vec3i
	{
		int x, y, z;
	};

	struct alignas( 16 ) Box
	{
		std::array<double, 3> xyz_min, xyz_max;
		bool contains( const vec3& pos ) const;
	};

	bool bboxes_overlap( const Box& a, const Box& b );
	void bbox_union( Box& rdi, const Box& a, const Box& b );

	// squared distance between a point inside the box, and the box surface
	double inner_point_box_squared_distance( const vec3& p, const Box& B );

	double tetra_signed_volume( const vec3& p1, const vec3& p2, const vec3& p3, const vec3& p4 );

	double point_triangle_squared_distance( const vec3& point, const vec3& V0, const vec3& V1, const vec3& V2, vec3& closest_point );

	inline double point_triangle_squared_distance( const vec3& point, const vec3& V0, const vec3& V1, const vec3& V2 )
	{
		vec3 cp;
		return point_triangle_squared_distance( point, V0, V1, V2, cp );
	}

	inline double geo_sqr( double d )
	{
		return d * d;
	}
}  // namespace GEO2