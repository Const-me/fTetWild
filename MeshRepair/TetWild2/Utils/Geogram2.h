#pragma once
#define CUSTOM_CODEZ 0
#if CUSTOM_CODEZ
#include "GeometricPrimitives.h"
#include "NearestSearch.h"
#include "orient3D.h"
#else
#include <geogram/basic/geometry.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/numerics/predicates.h>
#endif

namespace GEO2
{
	using coord_index_t = uint8_t;
	using index_t = uint32_t;

	constexpr uint32_t NO_FACET = ~( (uint32_t)0 );
#if !CUSTOM_CODEZ
	using GEO::vec3i;
	using GEO::vec3;
	using GEO::Box;
	using GEO::geo_sqr;
	using GEO::normalize;
	using GEO::distance;
	using GEO::dot;
	using GEO::length;
	using GEO::cross;
	using GEO::geo_sgn;
	using GEO::Geom::tetra_signed_volume;
	using GEO::PCK::orient_3d;

	inline double point_triangle_squared_distance( const vec3& point, const vec3& V0, const vec3& V1, const vec3& V2, vec3& closest_point )
	{
		double lambda0, lambda1, lambda2;
		return GEO::Geom::point_triangle_squared_distance( point, V0, V1, V2, closest_point, lambda0, lambda1, lambda2 );
	}

	inline double point_triangle_squared_distance( const vec3& point, const vec3& V0, const vec3& V1, const vec3& V2 )
	{
		vec3 closest_point;
		double lambda0, lambda1, lambda2;
		return GEO::Geom::point_triangle_squared_distance( point, V0, V1, V2, closest_point, lambda0, lambda1, lambda2 );
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
#endif
}  // namespace GEO2