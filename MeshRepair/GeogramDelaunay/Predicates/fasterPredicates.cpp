#include <geogram/basic/common.h>

// This makes sure the compiler will not optimize y = a*x+b
// with fused multiply-add, this would break the exact
// predicates.
#ifdef _MSC_VER
#pragma fp_contract( off )
#endif

#include <geogram/numerics/predicates.h>
#include <geogram/numerics/multi_precision.h>
#include <geogram/basic/assert.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/matrix.h>
#include <algorithm>

#include "RobustPredicates.h"
#define PCK_STAT( x )

namespace GEO
{
	namespace PCK
	{
		constexpr int FPG_UNCERTAIN_VALUE = 0;

		// The original header is only used for testing
#include <geogram/numerics/predicates/orient3d.h>

		Sign orient_3d_exact( const double* p0, const double* p1, const double* p2, const double* p3 )
		{
			PCK_STAT( cnt_orient3d_exact++ );

			const expansion& a11 = expansion_diff( p1[ 0 ], p0[ 0 ] );
			const expansion& a12 = expansion_diff( p1[ 1 ], p0[ 1 ] );
			const expansion& a13 = expansion_diff( p1[ 2 ], p0[ 2 ] );

			const expansion& a21 = expansion_diff( p2[ 0 ], p0[ 0 ] );
			const expansion& a22 = expansion_diff( p2[ 1 ], p0[ 1 ] );
			const expansion& a23 = expansion_diff( p2[ 2 ], p0[ 2 ] );

			const expansion& a31 = expansion_diff( p3[ 0 ], p0[ 0 ] );
			const expansion& a32 = expansion_diff( p3[ 1 ], p0[ 1 ] );
			const expansion& a33 = expansion_diff( p3[ 2 ], p0[ 2 ] );

			const expansion& Delta = expansion_det3x3( a11, a12, a13, a21, a22, a23, a31, a32, a33 );

			PCK_STAT( len_orient3d = std::max( len_orient3d, Delta.length() ) );

			return Delta.sign();
		}

		static Sign orient_3d_orig( const double* p0, const double* p1, const double* p2, const double* p3 )
		{
			PCK_STAT( cnt_orient3d_total++ );
			Sign result = Sign( orient_3d_filter( p0, p1, p2, p3 ) );
			if( result == 0 )
			{
				result = orient_3d_exact( p0, p1, p2, p3 );
			}
			return result;
		}

		inline Sign invertedSign( double r )
		{
			if( r < 0 )
				return Sign::POSITIVE;
			else if( r > 0 )
				return Sign::NEGATIVE;
			else
				return Sign::ZERO;
		}

		inline Sign orient3d_fast( const double* p0, const double* p1, const double* p2, const double* p3 )
		{
			const double r = RobustPredicates::orient3d( p0, p1, p2, p3 );
			return invertedSign( r );
		}

		Sign orient_3d( const double* p0, const double* p1, const double* p2, const double* p3 )
		{
#if 1
			return orient3d_fast( p0, p1, p2, p3 );
#else
			Sign fast = orient3d_fast( p0, p1, p2, p3 );
			Sign orig = orient_3d_orig( p0, p1, p2, p3 );
			if( orig == fast )
				return fast;
			__debugbreak();
			orient3d_fast( p0, p1, p2, p3 );
			orient_3d_orig( p0, p1, p2, p3 );
			return orig;
#endif
		}

		// The original header is only used for testing
#include <geogram/numerics/predicates/orient2d.h>

		static Sign orient_2d_exact( const double* p0, const double* p1, const double* p2 )
		{
			PCK_STAT( cnt_orient2d_exact++ );

			const expansion& a11 = expansion_diff( p1[ 0 ], p0[ 0 ] );
			const expansion& a12 = expansion_diff( p1[ 1 ], p0[ 1 ] );

			const expansion& a21 = expansion_diff( p2[ 0 ], p0[ 0 ] );
			const expansion& a22 = expansion_diff( p2[ 1 ], p0[ 1 ] );

			const expansion& Delta = expansion_det2x2( a11, a12, a21, a22 );

			PCK_STAT( len_orient2d = std::max( len_orient2d, Delta.length() ) );

			return Delta.sign();
		}

		static Sign orient_2d_orig( const double* p0, const double* p1, const double* p2 )
		{
			PCK_STAT( cnt_orient2d_total++ );
			Sign result = Sign( orient_2d_filter( p0, p1, p2 ) );
			if( result == 0 )
			{
				result = orient_2d_exact( p0, p1, p2 );
			}
			return result;
		}

		inline Sign makeSign( double r )
		{
			if( r < 0 )
				return Sign::NEGATIVE;
			else if( r > 0 )
				return Sign::POSITIVE;
			else
				return Sign::ZERO;
		}

		inline Sign orient2d_fast( const double* p0, const double* p1, const double* p2 )
		{
			const double r = RobustPredicates::orient2d( p0, p1, p2 );
			return makeSign( r );
		}

		Sign orient_2d( const double* p0, const double* p1, const double* p2 )
		{
#if 1
			return orient2d_fast( p0, p1, p2 );
#else
			Sign fast = orient2d_fast( p0, p1, p2 );
			Sign orig = orient_2d_orig( p0, p1, p2 );
			if( orig == fast )
				return fast;
			__debugbreak();
			return orig;
#endif
		}
	}  // namespace PCK
}  // namespace GEO