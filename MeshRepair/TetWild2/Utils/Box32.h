#pragma once
#include <array>
#include <immintrin.h>

namespace floatTetWild
{
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
	};
}  // namespace floatTetWild