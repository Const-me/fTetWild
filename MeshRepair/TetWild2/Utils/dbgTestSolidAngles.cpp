#include "stdafx.h"
#include "miscUtils.h"
#include "../includeEigen.h"
#include "AvxMath.h"

namespace
{
	constexpr double PI = 3.1415926535897932384626433832795;

	// Original version from IGL
	static double solidAngleOrig( const double* A, const double* B, const double* C, const double* P )
	{
		using SType = double;
		// Gather vectors to corners
		Eigen::Matrix<SType, 3, 3> v;
		// Don't use this since it will freak out for templates with != 3 size
		// v<< (A-P),(B-P),(C-P);
		for( int d = 0; d < 3; d++ )
		{
			v( 0, d ) = A[ d ] - P[ d ];
			v( 1, d ) = B[ d ] - P[ d ];
			v( 2, d ) = C[ d ] - P[ d ];
		}
		Eigen::Matrix<SType, 1, 3> vl = v.rowwise().norm();
		// printf("\n");
		//  Compute determinant
		SType detf = v( 0, 0 ) * v( 1, 1 ) * v( 2, 2 ) + v( 1, 0 ) * v( 2, 1 ) * v( 0, 2 ) + v( 2, 0 ) * v( 0, 1 ) * v( 1, 2 ) -
					 v( 2, 0 ) * v( 1, 1 ) * v( 0, 2 ) - v( 1, 0 ) * v( 0, 1 ) * v( 2, 2 ) - v( 0, 0 ) * v( 2, 1 ) * v( 1, 2 );
		// Compute pairwise dotproducts
		Eigen::Matrix<SType, 1, 3> dp;
		dp( 0 ) = v( 1, 0 ) * v( 2, 0 );
		dp( 0 ) += v( 1, 1 ) * v( 2, 1 );
		dp( 0 ) += v( 1, 2 ) * v( 2, 2 );
		dp( 1 ) = v( 2, 0 ) * v( 0, 0 );
		dp( 1 ) += v( 2, 1 ) * v( 0, 1 );
		dp( 1 ) += v( 2, 2 ) * v( 0, 2 );
		dp( 2 ) = v( 0, 0 ) * v( 1, 0 );
		dp( 2 ) += v( 0, 1 ) * v( 1, 1 );
		dp( 2 ) += v( 0, 2 ) * v( 1, 2 );
		// Compute winding number
		// Only divide by TWO_PI instead of 4*pi because there was a 2 out front
		return atan2( detf, vl( 0 ) * vl( 1 ) * vl( 2 ) + dp( 0 ) * vl( 0 ) + dp( 1 ) * vl( 1 ) + dp( 2 ) * vl( 2 ) ) / ( 2.0 * PI );
	}

	// Compute lengths of 3 vectors
	inline __m256d computeLengths( __m256d v0, __m256d v1, __m256d v2 )
	{
		// Reduce count of shuffles by transposing them.
		__m256d c0, c1, c2;
		AvxMath::transpose3x3( v0, v1, v2, c0, c1, c2 );

		c0 = _mm256_mul_pd( c0, c0 );
		c1 = _mm256_mul_pd( c1, c1 );
		c2 = _mm256_mul_pd( c2, c2 );

		__m256d len2 = _mm256_add_pd( _mm256_add_pd( c0, c1 ), c2 );
		return _mm256_sqrt_pd( len2 );
	}

	// Compute determinant of the 3x3 matrix
	inline double det3x3( __m256d v0, __m256d v1, __m256d v2 )
	{
		// https://en.wikipedia.org/wiki/Triple_product#Properties
		const __m256d cp = AvxMath::vector3Cross( v1, v2 );
		return AvxMath::vector3DotScalar( cp, v0 );
	}

	// Compute horizontal product of the 3D vector
	inline double vector3HorizontalProduct( __m256d v )
	{
		__m128d xy = AvxMath::low2( v );
		__m128d z = AvxMath::high2( v );
		xy = _mm_mul_sd( xy, _mm_unpackhi_pd( xy, xy ) );
		xy = _mm_mul_sd( xy, z );
		return _mm_cvtsd_f64( xy );
	}

	static double __declspec( noinline ) solidAngleAvx( const double* A, const double* B, const double* C, const double* P )
	{
		using namespace AvxMath;

		// Load all 4 inputs
		const __m256d a = loadDouble3( A );
		const __m256d b = loadDouble3( B );
		const __m256d c = loadDouble3( C );
		const __m256d p = loadDouble3( P );

		// Compute view directions
		const __m256d v0 = _mm256_sub_pd( a, p );
		const __m256d v1 = _mm256_sub_pd( b, p );
		const __m256d v2 = _mm256_sub_pd( c, p );

		const __m256d lengths = computeLengths( v0, v1, v2 );

		// Compute determinant
		const double detf = det3x3( v0, v1, v2 );

		// Compute pairwise dot products
		const double dp0 = vector3DotScalar( v1, v2 );
		const double dp1 = vector3DotScalar( v2, v0 );
		const double dp2 = vector3DotScalar( v0, v1 );
		const __m256d dp = _mm256_setr_pd( dp0, dp1, dp2, 0 );

		const double hp = vector3HorizontalProduct( lengths );
		const double dotdot = vector3DotScalar( lengths, dp );
		return std::atan2( detf, hp + dotdot ) / ( 2.0 * PI );
	}
}  // namespace

void dbgTestSolidAngles()
{
	std::array<double, 12> arr;
	srand( 0 );
	for( double& d : arr )
		d = (double)rand();

	const double orig = solidAngleOrig( &arr[ 0 ], &arr[ 3 ], &arr[ 6 ], &arr[ 9 ] );
	const double my = solidAngleAvx( &arr[ 0 ], &arr[ 3 ], &arr[ 6 ], &arr[ 9 ] );
	printf( "dbgTestSolidAngles: IGL %g, AVX %g\n", orig, my );
	__debugbreak();
}