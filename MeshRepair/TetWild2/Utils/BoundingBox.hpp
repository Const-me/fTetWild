#pragma once
#include "Geogram2.h"

// Utility class to compute bounding boxes of things with SSE
class BoundingBox
{
	__m128d xyMin, zMin, xyMax, zMax;

  public:
	// Initialize with numeric limit for doubles
	BoundingBox()
	{
		xyMin = _mm_set1_pd( DBL_MAX );
		zMin = xyMin;
		xyMax = _mm_sub_pd( _mm_setzero_pd(), xyMin );
		zMax = xyMax;
	}

	// Initialize min/max to the specified position
	BoundingBox( const double* rsi )
	{
		xyMin = _mm_loadu_pd( rsi );
		xyMax = xyMin;
		zMin = _mm_load_sd( rsi + 2 );
		zMax = zMin;
	}

	// Extend box with a point
	void extend( const double* rsi )
	{
		__m128d v = _mm_loadu_pd( rsi );
		xyMin = _mm_min_pd( xyMin, v );
		xyMax = _mm_max_pd( xyMax, v );

		v = _mm_load_sd( rsi + 2 );
		zMin = _mm_min_pd( zMin, v );
		zMax = _mm_max_pd( zMax, v );
	}

	// Compute diagonal of the box
	double diagonal() const
	{
		__m128d xy = _mm_sub_pd( xyMax, xyMin );
		__m128d z = _mm_sub_pd( zMax, zMin );

		xy = _mm_dp_pd( xy, xy, 0b00110001 );
		z = _mm_mul_sd( z, z );

		__m128d res = _mm_add_sd( xy, z );
		res = _mm_sqrt_sd( res, res );
		return _mm_cvtsd_f64( res );
	}

	// Store the box as FP64 values in memory
	void store( double* pMin, double* pMax ) const
	{
		_mm_storeu_pd( pMin, xyMin );
		_mm_store_sd( pMin + 2, zMin );

		_mm_storeu_pd( pMax, xyMax );
		_mm_store_sd( pMax + 2, zMax );
	}

	void store( GEO2::Box& rdi ) const
	{
		// TODO [low]: rework into 3 aligned stores, using unpacklo and blend to combine vectors
		store( &rdi.xyz_min[ 0 ], &rdi.xyz_max[ 0 ] );
	}
};