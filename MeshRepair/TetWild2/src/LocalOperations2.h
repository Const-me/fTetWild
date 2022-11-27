#pragma once
#include <array>
#include "Types.hpp"

namespace floatTetWild
{
	// This version is only slightly faster than the original, but the output is bitwise identical
	double AMIPS_energy_aux_v2( const std::array<double, 12>& arr );

	// Minor math optimizations, breaks binary equality with the original
	double AMIPS_energy_aux_v3( const std::array<double, 12>& arr );

	// Vectorized into AVX
	double AMIPS_energy_aux_v4( const std::array<double, 12>& arr );

	// Optimized the energy formula for better numerical accuracy
	double AMIPS_energy_aux_v5( const std::array<double, 12>& arr );

	// Vectorized version which is still bitwise identical to the original one
	void AMIPS_hessian_v2( const std::array<Scalar, 12>& T, Matrix3& result_0 );

	// Vectorized version which no longer bitwise identical, due to different summation orders and minor algebra optimizations
	void AMIPS_hessian_v3( const std::array<Scalar, 12>& T, Matrix3& result_0 );

	// Vectorized version which translates vertices to barycenter, for better numerical accuracy
	void AMIPS_hessian_v4( const std::array<Scalar, 12>& T, Matrix3& result_0 );

	void AMIPS_jacobian_v2( const std::array<Scalar, 12>& T, Vector3& result_0 );
}