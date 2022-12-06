// This file is part of TriWild, a software for generating linear/curved triangle meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "mpir/gmp.h"
#include <iostream>

namespace triwild
{
	// An infinite-precision rational number implemented in that mpir.dll third-party library.
	class Rational
	{
		// Based on the source code in clear.c and memory.c, the default implementation calls free() in release builds
		// This means when all these pointers are 0, the mpq_clear() function gonna do nothing, safely so.
		void setZero()
		{
			static_assert( sizeof( __mpz_struct ) == 16 );
			const __m128i z = _mm_setzero_si128();
			_mm_storeu_si128( (__m128i*)&value->_mp_num, z );
			_mm_storeu_si128( (__m128i*)&value->_mp_den, z );
		}

		mpq_t value;

	  public:
		// Remove common factors of the denominator and numerator
		void canonicalize()
		{
			mpq_canonicalize( value );
		}

		// Return +1 if this value is a positive number, -1 if this value is negative, and zero if this value is 0
		int getSign() const
		{
			return mpq_sgn( value );
		}

		Rational()
		{
			// mpq_init - Make a new rational number with value 0/1, in init.c
			mpq_init( value );
		}

		Rational( double d )
		{
			mpq_init( value );
			mpq_set_d( value, d );
			// canonicalize();
		}

		// Make a rational number numerator/denominator
		Rational( int numerator, int denominator )
		{
			assert( 0 != denominator );
			mpq_init( value );
			mpq_set_si( value, numerator, denominator );
		}

		// Set this value to numerator / denominator
		void rational( int numerator, int denominator )
		{
			assert( 0 != denominator );
			mpq_set_si( value, numerator, denominator );
		}

		Rational( const Rational& other )
		{
			mpq_init( value );
			mpq_set( value, other.value );
		}

		Rational( Rational&& other ) noexcept
		{
			value[ 0 ] = other.value[ 0 ];
			other.setZero();
		}

		~Rational()
		{
			mpq_clear( value );
		}

		friend Rational operator-( const Rational& x )
		{
			Rational r_out;
			mpq_neg( r_out.value, x.value );
			return r_out;
		}

		friend Rational operator+( const Rational& x, const Rational& y )
		{
			Rational r_out;
			mpq_add( r_out.value, x.value, y.value );
			return r_out;
		}

		friend Rational operator-( const Rational& x, const Rational& y )
		{
			Rational r_out;
			mpq_sub( r_out.value, x.value, y.value );
			return r_out;
		}
		friend Rational operator*( const Rational& x, const Rational& y )
		{
			Rational r_out;
			mpq_mul( r_out.value, x.value, y.value );
			return r_out;
		}

		void operator+=( const Rational& that )
		{
			// See aors.c, the function supports accumulation
			mpq_add( value, value, that.value );
		}
		void operator-=( const Rational& that )
		{
			// See aors.c, the function supports accumulation
			mpq_sub( value, value, that.value );
		}
		void operator*=( const Rational& that )
		{
			// See mul.c, they support accumulation there
			mpq_mul( value, value, that.value );
		}
		void operator/=( const Rational& that )
		{
			// See div.c, they support accumulation there
			mpq_div( value, value, that.value );
		}
		// Set this value to a + b
		void add( const Rational& a, const Rational& b )
		{
			mpq_add( value, a.value, b.value );
		}
		// Set this value to a - b
		void sub( const Rational& a, const Rational& b )
		{
			mpq_sub( value, a.value, b.value );
		}
		// Set this value to a * b
		void mul( const Rational& a, const Rational& b )
		{
			mpq_mul( value, a.value, b.value );
		}
		void negate()
		{
			mpq_neg( value, value );
		}

		friend Rational operator/( const Rational& x, const Rational& y )
		{
			Rational r_out;
			mpq_div( r_out.value, x.value, y.value );
			return r_out;
		}

		static inline Rational zero()
		{
			return Rational {};
		}
		static inline Rational one()
		{
			return Rational { 1, 1 };
		}

		friend Rational pow( const Rational& x, int p )
		{
			if( 0 != p )
			{
				Rational r_out = x;
				for( int i = 1; i < std::abs( p ); i++ )
					r_out *= x;

				if( p < 0 )
				{
					// Based on the source code in inv.c, that function supports passing same pointers to both arguments
					mpq_inv( r_out.value, r_out.value );
				}
				return r_out;
			}
			return one();
		}

		void operator=( const Rational& x )
		{
			if( this == &x )
				return;
			mpq_set( value, x.value );
		}

		void operator=( Rational&& x ) noexcept
		{
			if( this == &x )
				return;
			mpq_swap( value, x.value );
		}

		void operator=( double x )
		{
			mpq_set_d( value, x );
			// canonicalize();
		}

		//> < ==
		friend bool operator<( const Rational& r, const Rational& r1 )
		{
			return mpq_cmp( r.value, r1.value ) < 0;
		}

		friend bool operator>( const Rational& r, const Rational& r1 )
		{
			return mpq_cmp( r.value, r1.value ) > 0;
		}

		friend bool operator<=( const Rational& r, const Rational& r1 )
		{
			return mpq_cmp( r.value, r1.value ) <= 0;
		}

		friend bool operator>=( const Rational& r, const Rational& r1 )
		{
			return mpq_cmp( r.value, r1.value ) >= 0;
		}

		friend bool operator==( const Rational& r, const Rational& r1 )
		{
			return mpq_equal( r.value, r1.value );
		}

		friend bool operator!=( const Rational& r, const Rational& r1 )
		{
			return !mpq_equal( r.value, r1.value );
		}

		// Convert to FP64, using rounding towards zero
		double asDouble() const
		{
			return mpq_get_d( value );
		}

		//<<
		friend std::ostream& operator<<( std::ostream& os, const Rational& r )
		{
			os << mpq_get_d( r.value );
			return os;
		}
	};
}  // namespace triwild