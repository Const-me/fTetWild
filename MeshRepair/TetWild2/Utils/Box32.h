#pragma once
#include <array>
#include <immintrin.h>

namespace floatTetWild
{
	// Downcast FP64 vector to FP32, rounding towards negative infinity
	__m128 downcastFloor( __m256d pos );

	// Downcast FP64 vector to FP32, rounding towards positive infinity
	__m128 downcastCeil( __m256d pos );

	// Axis-aligned 3D bounding box, with FP32 precision
	struct alignas( 8 ) Box32
	{
		std::array<float, 3> boxMin, boxMax;

		// Initialize box for the triangle, using proper rounding from FP64 to FP32
		void computeTriangle( __m256d a, __m256d b, __m256d c );

		// Compute union of 2 bounding boxes
		void computeUnion( const Box32& a, const Box32& b );

		// Compute squared distance between a point and this Box.
		// When the point is inside, the result is negative, equal to the distance^2 to the closest face of the box
		// When the point is outside, the result is positive, distance^2 to the box
		// The result is broadcasted into all 4 lanes of the output vector
		__m128 pointBoxSignedSquaredDistance( __m128 pos ) const;

		// Same as above, for 2 boxes and 2 input points in parallel
		static __m256 pointBoxSignedSquaredDistanceX2( const Box32& b1, const Box32& b2, __m256 pos );

		// Set this box to 0.0 in all 6 scalars
		void setZero();

		static __m256 createBoxVector( __m256d boxMin, __m256d boxMax );

		// True when this box intersects other box, passed in the input registers
		bool intersects( __m256 boxVector ) const;
	};
}  // namespace floatTetWild