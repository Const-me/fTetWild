#pragma once
#include <geogram/basic/geometry.h>
#include <geogram/numerics/predicates.h>
#ifdef __AVX__
#include "AvxMath.h"
#else
#error This code requires at least AVX1
#endif

namespace GEO2
{
	using coord_index_t = uint8_t;
	using index_t = uint32_t;

	constexpr uint32_t NO_FACET = ~( (uint32_t)0 );
	using GEO::Box;
	using GEO::cross;
	using GEO::distance;
	using GEO::dot;
	using GEO::geo_sgn;
	using GEO::geo_sqr;
	using GEO::length;
	using GEO::normalize;
	using GEO::vec3;
	using GEO::vec3i;
	using GEO::Geom::tetra_signed_volume;
	using GEO::PCK::orient_3d;

	double point_triangle_squared_distance( __m256d point, const vec3& V0, const vec3& V1, const vec3& V2, vec3* closest_point = nullptr );

	inline double point_triangle_squared_distance( const vec3& point, const vec3& V0, const vec3& V1, const vec3& V2, vec3* closest_point = nullptr )
	{
		__m256d p = AvxMath::loadDouble3( &point.x );
		return point_triangle_squared_distance( p, V0, V1, V2, closest_point );
	}

	inline double inner_point_box_squared_distance( const vec3& p, const Box& B )
	{
		geo_debug_assert( B.contains( p ) );
		double result = geo_sqr( p[ 0 ] - B.xyz_min[ 0 ] );
		result = std::min( result, geo_sqr( p[ 0 ] - B.xyz_max[ 0 ] ) );
		for( coord_index_t c = 1; c < 3; ++c )
		{
			result = std::min( result, geo_sqr( p[ c ] - B.xyz_min[ c ] ) );
			result = std::min( result, geo_sqr( p[ c ] - B.xyz_max[ c ] ) );
		}
		return result;
	}

	void dbgRunSomeTests();

	inline double distance2( __m256d a, const vec3& b )
	{
		using namespace AvxMath;
		__m256d bv = loadDouble3( &b.x );
		__m256d dist = _mm256_sub_pd( bv, a );
		return vector3DotScalar( dist, dist );
	}
}  // namespace GEO2