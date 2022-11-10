#pragma once

namespace RobustPredicates
{
	void exactinit();

	// A positive value if the points pa, pb, and pc occur in counterclockwise order
	// A negative value if they occur in clockwise order
	// Zero if they are collinear.
	// The result is also a rough approximation of twice the signed area of the triangle defined by the three points.
	double orient2d( const double* pa, const double* pb, const double* pc );

	// A positive value if the point pd lies below the plane passing through pa, pb, and pc;
	// "below" is defined so that pa, pb, and pc appear in counterclockwise order when viewed from above the plane.
	// Returns a negative value if pd lies above the plane.
	// Returns zero if the points are  coplanar.
	// The result is also a rough approximation of six times the signed volume of the tetrahedron defined by the four points.
	double orient3d( const double* pa, const double* pb, const double* pc, const double* pd );

	// A positive value if the point pd lies inside the circle passing through pa, pb, and pc
	// A negative value if it lies outside
	// Zero if the four points are cocircular.
	// The points pa, pb, and pc must be in counterclockwise order, or the sign of the result will be reversed.
	double incircle( const double* pa, const double* pb, const double* pc, const double* pd );

	// A positive value if the point pe lies inside the sphere passing through pa, pb, pc, and pd;
	// A negative value  if it lies outside;
	// Zero if the five points are cospherical.
	// The points pa, pb, pc, and pd must be ordered so that they have a positive orientation as defined by orient3d(), or the sign of the result will be reversed.
	double insphere( const double* pa, const double* pb, const double* pc, const double* pd, const double* pe );
}  // namespace RodustPredicates