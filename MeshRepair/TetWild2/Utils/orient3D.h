#pragma once
#include "GeometricPrimitives.h"

namespace GEO2
{
	/**
	 * \brief Computes the orientation predicate in 3d.
	 * \details Computes the sign of the signed volume of
	 *  the tetrahedron p0, p1, p2, p3.
	 * \param[in] p0 , p1 , p2 , p3 vertices of the tetrahedron
	 * \retval POSITIVE if the tetrahedron is oriented positively
	 * \retval ZERO if the tetrahedron is flat
	 * \retval NEGATIVE if the tetrahedron is oriented negatively
	 * \todo check whether orientation is inverted as compared to
	 *   Shewchuk's version.
	 */
	Sign orient_3d( const double* p0, const double* p1, const double* p2, const double* p3 );

	inline Sign orient_3d( const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3 )
	{
		return orient_3d( &p0.x, &p1.x, &p2.x, &p3.x );
	}
}