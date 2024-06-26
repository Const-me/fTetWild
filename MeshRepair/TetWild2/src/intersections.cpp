// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "stdafx.h"
#include "intersections.h"
#include "LocalOperations.h"
#include "../Utils/triangleIntersection.h"

bool floatTetWild::seg_line_intersection_2d( const std::array<Vector2, 2>& seg, const std::array<Vector2, 2>& line, Scalar& t_seg )
{
	// assumptions:
	// segs are not degenerate
	// not coplanar

	const Scalar& x1 = seg[ 0 ][ 0 ];
	const Scalar& y1 = seg[ 0 ][ 1 ];
	const Scalar& x2 = seg[ 1 ][ 0 ];
	const Scalar& y2 = seg[ 1 ][ 1 ];

	const Scalar& x3 = line[ 0 ][ 0 ];
	const Scalar& y3 = line[ 0 ][ 1 ];
	const Scalar& x4 = line[ 1 ][ 0 ];
	const Scalar& y4 = line[ 1 ][ 1 ];

	Scalar n1 = ( y3 - y4 ) * ( x1 - x3 ) + ( x4 - x3 ) * ( y1 - y3 );
	Scalar d1 = ( x4 - x3 ) * ( y1 - y2 ) - ( x1 - x2 ) * ( y4 - y3 );
	if( d1 == 0 )
		return false;
	t_seg = n1 / d1;
	//    Scalar n2 = (y1 - y2) * (x1 - x3) + (x2 - x1) * (y1 - y3);
	Scalar d2 = ( x4 - x3 ) * ( y1 - y2 ) - ( x1 - x2 ) * ( y4 - y3 );
	if( d2 == 0 )
	{
		//        cout<<"d2==0"<<endl;
		return false;
	}

	if( t_seg < 0 || t_seg > 1 )
	{
		//        cout<<"t_seg = "<<t_seg<<endl;
		return false;
	}

	return true;
}

bool floatTetWild::seg_seg_intersection_2d( const std::array<Vector2, 2>& seg1, const std::array<Vector2, 2>& seg2, Scalar& t2 )
{
	// assumptions:
	// segs are not degenerate
	// not coplanar

	const Scalar& x1 = seg1[ 0 ][ 0 ];
	const Scalar& y1 = seg1[ 0 ][ 1 ];
	const Scalar& x2 = seg1[ 1 ][ 0 ];
	const Scalar& y2 = seg1[ 1 ][ 1 ];

	const Scalar& x3 = seg2[ 0 ][ 0 ];
	const Scalar& y3 = seg2[ 0 ][ 1 ];
	const Scalar& x4 = seg2[ 1 ][ 0 ];
	const Scalar& y4 = seg2[ 1 ][ 1 ];

	Scalar n1 = ( y3 - y4 ) * ( x1 - x3 ) + ( x4 - x3 ) * ( y1 - y3 );
	Scalar d1 = ( x4 - x3 ) * ( y1 - y2 ) - ( x1 - x2 ) * ( y4 - y3 );
	if( d1 == 0 )
		return false;
	Scalar t1 = n1 / d1;
	Scalar n2 = ( y1 - y2 ) * ( x1 - x3 ) + ( x2 - x1 ) * ( y1 - y3 );
	Scalar d2 = ( x4 - x3 ) * ( y1 - y2 ) - ( x1 - x2 ) * ( y4 - y3 );
	if( d2 == 0 )
		return false;
	t2 = n2 / d2;

	if( t1 < 0 || t1 > 1 || t2 < 0 || t2 > 1 )
		return false;

	return true;
}

floatTetWild::Scalar floatTetWild::seg_seg_squared_dist_3d( const std::array<Vector3, 2>& s1, const std::array<Vector3, 2>& s2 )
{
	Vector3 w0 = s1[ 0 ] - s2[ 0 ];
	Vector3 u = ( s1[ 1 ] - s1[ 0 ] ).normalized();
	Vector3 v = ( s2[ 1 ] - s2[ 0 ] ).normalized();
	Scalar a = u.dot( u );
	Scalar b = u.dot( v );
	Scalar c = v.dot( v );
	Scalar d = u.dot( w0 );
	Scalar e = v.dot( w0 );

	Scalar dd = a * c - b * b;
	Scalar t1, t2;
	if( dd == 0 )
	{
		t1 = 0;
		t2 = d / b;
	}
	else
	{
		t1 = ( b * e - c * d ) / dd;
		t2 = ( a * e - b * d ) / dd;
	}

	return ( ( ( 1 - t1 ) * s1[ 0 ] + t1 * s1[ 1 ] ) - ( ( 1 - t2 ) * s2[ 0 ] + t2 * s2[ 1 ] ) ).squaredNorm();
}

floatTetWild::Scalar floatTetWild::p_line_squared_dist_3d( const Vector3& v, const Vector3& a, const Vector3& b )
{
	return ( ( b - a ).cross( a - v ) ).squaredNorm() / ( b - a ).squaredNorm();
}

floatTetWild::Scalar floatTetWild::p_seg_squared_dist_3d( const Vector3& v, const Vector3& a, const Vector3& b )
{
	Vector3 av = v - a;
	Vector3 ab = b - a;
	if( av.dot( ab ) < 0 )
		return av.squaredNorm();
	Vector3 bv = v - b;
	if( bv.dot( -ab ) < 0 )
		return bv.squaredNorm();

	return ( ab.cross( -av ) ).squaredNorm() / ab.squaredNorm();
}

bool floatTetWild::seg_plane_intersection( const Vector3& p1, const Vector3& p2, const Vector3& a, const Vector3& n, Vector3& p, Scalar& d1 )
{
	Vector3 u = p2 - p1;
	Vector3 w = p1 - a;

	Scalar D = n.dot( u );
	d1 = -n.dot( w );

	//    if (fabs(D) <= SCALAR_ZERO) {// segment is parallel to plane
	Scalar t;
	if( fabs( D ) == 0 )
	{				   // segment is parallel to plane
		if( d1 == 0 )  // segment lies in plane
			t = INT_MAX;
		else
			t = INT_MIN;
		//            return 2;
		//        else // no intersection
		//            return 0;
		return false;
	}

	t = d1 / D;
	if( t <= 0 || t >= 1 )
	{
		return false;
	}

	p = p1 + t * u;
	return true;
}

floatTetWild::eCutResult floatTetWild::is_tri_tri_cutted(
  const std::array<Vector3, 3>& f_tri, const std::array<Vector3, 3>& f_tet, const std::array<eOrientation, 3>& oris_tri )
{
	int cnt_pos = 0;
	int cnt_neg = 0;
	int cnt_zero = 0;
	for( int j = 0; j < 3; j++ )
	{
		if( oris_tri[ j ] == eOrientation::Zero )
			cnt_zero++;
		else if( oris_tri[ j ] == eOrientation::Positive )
			cnt_pos++;
		else
			cnt_neg++;
	}
	if( cnt_pos == 3 || cnt_neg == 3 )
		return eCutResult::Empty;

	if( cnt_zero == 3 )
	{
		int t = get_t( f_tri[ 0 ], f_tri[ 1 ], f_tri[ 2 ] );
		if( is_tri_tri_cutted_2d( { { to_2d( f_tri[ 0 ], t ), to_2d( f_tri[ 1 ], t ), to_2d( f_tri[ 2 ], t ) } },
			  { { to_2d( f_tet[ 0 ], t ), to_2d( f_tet[ 1 ], t ), to_2d( f_tet[ 2 ], t ) } } ) )
			return eCutResult::Coplanar;
		return eCutResult::Empty;
	}

	if( cnt_zero == 2 )
	{
		int t = get_t( f_tri[ 0 ], f_tri[ 1 ], f_tri[ 2 ] );
		for( int j = 0; j < 3; j++ )
		{
			if( oris_tri[ j ] == oris_tri[ ( j + 1 ) % 3 ] && is_seg_tri_cutted_2d( { { to_2d( f_tet[ j ], t ), to_2d( f_tet[ ( j + 1 ) % 3 ], t ) } },
																{ { to_2d( f_tri[ 0 ], t ), to_2d( f_tri[ 1 ], t ), to_2d( f_tri[ 2 ], t ) } } ) )
			{
				if( j == 0 )
					return eCutResult::Edge0;
				if( j == 1 )
					return eCutResult::Edge1;
				if( j == 2 )
					return eCutResult::Edge2;
			}
		}
		return eCutResult::Empty;
	}

	if( cnt_pos == 0 || cnt_neg == 0 )
		return eCutResult::Empty;

	std::array<eOrientation, 3> oris_tet;
	int cnt_pos1 = 0;
	int cnt_neg1 = 0;
	int cnt_zero1 = 0;
	for( int j = 0; j < 3; j++ )
	{
		oris_tet[ j ] = Predicates::orient_3d( f_tet[ 0 ], f_tet[ 1 ], f_tet[ 2 ], f_tri[ j ] );
		if( oris_tet[ j ] == eOrientation::Zero )
			cnt_zero1++;
		else if( oris_tet[ j ] == eOrientation::Positive )
			cnt_pos1++;
		else
			cnt_neg1++;
	}
	if( cnt_pos1 == 0 || cnt_neg1 == 0 )
		return eCutResult::Empty;

	auto is_f1_cut_f2 = []( const std::array<Vector3, 3>& f1, const std::array<Vector3, 3>& f2, const std::array<eOrientation, 3>& oris ) {	 // check f2 plane
		Vector3 n = ( ( f2[ 1 ] - f2[ 2 ] ).cross( f2[ 0 ] - f2[ 2 ] ) ).normalized();
		int t = get_t( f2[ 0 ], f2[ 1 ], f2[ 2 ] );
		for( int j = 0; j < 3; j++ )
		{
			if( oris[ j ] == eOrientation::Negative && oris[ ( j + 1 ) % 3 ] == eOrientation::Positive ||
				oris[ ( j + 1 ) % 3 ] == eOrientation::Negative && oris[ j ] == eOrientation::Positive )
			{
				Vector3 p;
				Scalar _;
				if( seg_plane_intersection( f1[ j ], f1[ ( j + 1 ) % 3 ], n, f2[ 0 ], p, _ ) &&
					is_p_inside_tri_2d( to_2d( p, t ), { { to_2d( f2[ 0 ], t ), to_2d( f2[ 1 ], t ), to_2d( f2[ 2 ], t ) } } ) )
					return true;
			}
		}
		return false;
	};

	// check tri plane
	if( is_f1_cut_f2( f_tet, f_tri, oris_tri ) )
		return eCutResult::Face;

	// check tet plane
	if( is_f1_cut_f2( f_tri, f_tet, oris_tet ) )
		return eCutResult::Face;

	return eCutResult::Empty;
}

bool floatTetWild::is_tri_tri_cutted_2d( const std::array<Vector2, 3>& vs_tet, const std::array<Vector2, 3>& vs_tri )
{
	std::array<eOrientation, 9> tri_tet;
	int cnt_pos0 = 0;
	int cnt_neg0 = 0;
	int cnt_pos1 = 0;
	int cnt_neg1 = 0;
	int cnt_pos2 = 0;
	int cnt_neg2 = 0;
	for( int i = 0; i < 3; i++ )
	{
		tri_tet[ i * 3 ] = Predicates::orient_2d( vs_tri[ i ], vs_tri[ ( i + 1 ) % 3 ], vs_tet[ 0 ] );
		if( tri_tet[ i * 3 ] == eOrientation::Positive )
			cnt_pos0++;
		else if( tri_tet[ i * 3 ] == eOrientation::Negative )
			cnt_neg0++;

		tri_tet[ i * 3 + 1 ] = Predicates::orient_2d( vs_tri[ i ], vs_tri[ ( i + 1 ) % 3 ], vs_tet[ 1 ] );
		if( tri_tet[ i * 3 + 1 ] == eOrientation::Positive )
			cnt_pos1++;
		else if( tri_tet[ i * 3 + 1 ] == eOrientation::Negative )
			cnt_neg1++;

		tri_tet[ i * 3 + 2 ] = Predicates::orient_2d( vs_tri[ i ], vs_tri[ ( i + 1 ) % 3 ], vs_tet[ 2 ] );
		if( tri_tet[ i * 3 + 2 ] == eOrientation::Positive )
			cnt_pos2++;
		else if( tri_tet[ i * 3 + 2 ] == eOrientation::Negative )
			cnt_neg2++;
	}
	if( cnt_neg0 == 3 || cnt_pos0 == 3 || cnt_neg1 == 3 || cnt_pos1 == 3 || cnt_neg2 == 3 || cnt_pos2 == 3 )
		return true;  // one of the vertices is strictly contained inside the other tri
	if( std::find( tri_tet.begin(), tri_tet.end(), eOrientation::Negative ) == tri_tet.end() ||
		std::find( tri_tet.begin(), tri_tet.end(), eOrientation::Positive ) == tri_tet.end() )
		return true;  // tet face is strictly contained by tri face

	std::array<eOrientation, 9> tet_tri;
	cnt_neg0 = cnt_pos0 = cnt_neg1 = cnt_pos1 = cnt_neg2 = cnt_pos2 = 0;
	for( int i = 0; i < 3; i++ )
	{
		tet_tri[ i * 3 ] = Predicates::orient_2d( vs_tet[ i ], vs_tet[ ( i + 1 ) % 3 ], vs_tri[ 0 ] );
		if( tet_tri[ i * 3 ] == eOrientation::Positive )
			cnt_pos0++;
		else if( tet_tri[ i * 3 ] == eOrientation::Negative )
			cnt_neg0++;

		tet_tri[ i * 3 + 1 ] = Predicates::orient_2d( vs_tet[ i ], vs_tet[ ( i + 1 ) % 3 ], vs_tri[ 1 ] );
		if( tet_tri[ i * 3 + 1 ] == eOrientation::Positive )
			cnt_pos1++;
		else if( tet_tri[ i * 3 + 1 ] == eOrientation::Negative )
			cnt_neg1++;

		tet_tri[ i * 3 + 2 ] = Predicates::orient_2d( vs_tet[ i ], vs_tet[ ( i + 1 ) % 3 ], vs_tri[ 2 ] );
		if( tet_tri[ i * 3 + 2 ] == eOrientation::Positive )
			cnt_pos2++;
		else if( tet_tri[ i * 3 + 2 ] == eOrientation::Negative )
			cnt_neg2++;
	}
	if( cnt_neg0 == 3 || cnt_pos0 == 3 || cnt_neg1 == 3 || cnt_pos1 == 3 || cnt_neg2 == 3 || cnt_pos2 == 3 )
		return true;
	if( std::find( tet_tri.begin(), tet_tri.end(), eOrientation::Negative ) == tet_tri.end() ||
		std::find( tet_tri.begin(), tet_tri.end(), eOrientation::Positive ) == tet_tri.end() )
		return true;  // tri face is contained by tet face

	for( int tri_e_id = 0; tri_e_id < 3; tri_e_id++ )
	{
		for( int tet_e_id = 0; tet_e_id < 3; tet_e_id++ )
		{
			if( is_crossing( tri_tet[ tri_e_id * 3 + tet_e_id ], tri_tet[ tri_e_id * 3 + ( tet_e_id + 1 ) % 3 ] ) &&
				is_crossing( tet_tri[ tet_e_id * 3 + tri_e_id ], tet_tri[ tet_e_id * 3 + ( tri_e_id + 1 ) % 3 ] ) )
				return true;
		}
	}

	return false;
}

bool floatTetWild::is_seg_tri_cutted_2d( const std::array<Vector2, 2>& seg, const std::array<Vector2, 3>& tri )
{
	std::vector<Vector2> ps;
	for( int i = 0; i < 3; i++ )
	{
		Scalar t2;
		if( seg_seg_intersection_2d( seg, { { tri[ i ], tri[ ( i + 1 ) % 3 ] } }, t2 ) )
		{
			ps.push_back( t2 * tri[ i ] + ( 1 - t2 ) * tri[ ( i + 1 ) % 3 ] );	// todo double check
		}
	}
	if( ps.size() < 2 )
		return false;
	for( int i = 0; i < ps.size(); i++ )
	{
		Vector2 v = ps[ i ] - ps[ ( i + 1 ) % ps.size() ];
		if( v[ 0 ] < SCALAR_ZERO && v[ 1 ] < SCALAR_ZERO )
		{
			ps.erase( ps.begin() + i );
			i--;
		}
	}
	assert( ps.size() < 3 );
	if( ps.size() < 2 )
		return false;

	return true;

	//////////

	std::array<eOrientation, 3> tri_seg;
	for( int i = 0; i < 3; i++ )
	{
		tri_seg[ i ] = Predicates::orient_2d( seg[ 0 ], seg[ 1 ], tri[ i ] );
	}
	if( std::find( tri_seg.begin(), tri_seg.end(), eOrientation::Negative ) == tri_seg.end() ||
		std::find( tri_seg.begin(), tri_seg.end(), eOrientation::Positive ) == tri_seg.end() )
		return false;  // tri on the same side of seg

	std::array<eOrientation, 6> seg_tri;
	int cnt_pos0 = 0;
	int cnt_neg0 = 0;
	int cnt_pos1 = 0;
	int cnt_neg1 = 0;
	for( int i = 0; i < 3; i++ )
	{
		seg_tri[ i * 2 ] = Predicates::orient_2d( tri[ i ], tri[ ( i + 1 ) % 3 ], seg[ 0 ] );
		if( seg_tri[ i * 2 ] == eOrientation::Positive )
			cnt_pos0++;
		else if( seg_tri[ i * 2 ] == eOrientation::Negative )
			cnt_neg0++;
		seg_tri[ i * 2 + 1 ] = Predicates::orient_2d( tri[ i ], tri[ ( i + 1 ) % 3 ], seg[ 1 ] );
		if( seg_tri[ i * 2 + 1 ] == eOrientation::Positive )
			cnt_pos1++;
		else if( seg_tri[ i * 2 + 1 ] == eOrientation::Negative )
			cnt_neg1++;
	}
	if( cnt_neg0 == 3 || cnt_pos0 == 3 || cnt_neg1 == 3 || cnt_pos1 == 3 )	// one of the endpoints is contained inside tri
		return true;

	for( int tri_e_id = 0; tri_e_id < 3; tri_e_id++ )
	{
		if( is_crossing( seg_tri[ tri_e_id ], seg_tri[ ( tri_e_id + 1 ) % 3 ] ) && is_crossing( tri_seg[ tri_e_id ], tri_seg[ ( tri_e_id + 1 ) % 3 ] ) )
			return true;
	}

	return false;
}

bool floatTetWild::is_p_inside_tri_2d( const Vector2& p, const std::array<Vector2, 3>& tri )
{
	int cnt_pos = 0;
	int cnt_neg = 0;

	for( int i = 0; i < 3; i++ )
	{
		eOrientation ori = Predicates::orient_2d( tri[ i ], tri[ ( i + 1 ) % 3 ], p );
		if( ori == eOrientation::Positive )
			cnt_pos++;
		else if( ori == eOrientation::Negative )
			cnt_neg++;
	}
	//    if(cnt_neg==0 || cnt_pos==0)
	if( cnt_neg == 3 || cnt_pos == 3 )	// strict inside
		return true;
	return false;
}

int floatTetWild::get_t( const Vector3& p0, const Vector3& p1, const Vector3& p2 )
{
	using namespace AvxMath;
	// Load all 3 vectors
	const __m256d v0 = loadDouble3( p0.data() );
	const __m256d v1 = loadDouble3( p1.data() );
	const __m256d v2 = loadDouble3( p2.data() );

	// Compute the cross product
	__m256d n = vector3Cross( _mm256_sub_pd( v1, v2 ), _mm256_sub_pd( v0, v2 ) );
	// Absolute value
	n = vectorAbs( n );
	// Compute maximum of XYZ lanes
	const __m256d max4 = vector3BroadcastMaximum( n );
	// Return index of the first maximum lane
	return firstEqualLaneIndex( n, max4 );
}

int floatTetWild::get_t( const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, int idx )
{
	const Vector3i& tri = ib[ idx ];
	return get_t( vb[ tri[ 0 ] ], vb[ tri[ 1 ] ], vb[ tri[ 2 ] ] );
}

floatTetWild::Vector2 floatTetWild::to_2d( const Vector3& p, int t )
{
	return Vector2( p[ ( t + 1 ) % 3 ], p[ ( t + 2 ) % 3 ] );
}

floatTetWild::Vector2 floatTetWild::to_2d( const Vector3& p, const Vector3& n, const Vector3& pp, int t )
{
	Scalar dist = n.dot( p - pp );
	Vector3 proj_p = p - dist * n;
	return Vector2( proj_p[ ( t + 1 ) % 3 ], proj_p[ ( t + 2 ) % 3 ] );
}

bool floatTetWild::is_crossing( eOrientation s1, eOrientation s2 )
{
	if( s1 == eOrientation::Positive && s2 == eOrientation::Negative || s2 == eOrientation::Positive && s1 == eOrientation::Negative )
		return true;
	return false;
}

floatTetWild::eCutResult floatTetWild::is_tri_tri_cutted(
  const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& q1, const Vector3& q2, const Vector3& q3 )
{
	std::array<Scalar, 3> p_1 = { { 0, 0, 0 } }, q_1 = { { 0, 0, 0 } }, r_1 = { { 0, 0, 0 } };
	std::array<Scalar, 3> p_2 = { { 0, 0, 0 } }, q_2 = { { 0, 0, 0 } }, r_2 = { { 0, 0, 0 } };
	bool coplanar = false;
	std::array<Scalar, 3> s = { { 0, 0, 0 } }, t = { { 0, 0, 0 } };
	for( int j = 0; j < 3; j++ )
	{
		p_1[ j ] = p1[ j ];
		q_1[ j ] = p2[ j ];
		r_1[ j ] = p3[ j ];
		p_2[ j ] = q1[ j ];
		q_2[ j ] = q2[ j ];
		r_2[ j ] = q3[ j ];
	}

	if( !triangleIntersectionTest( &p_1[ 0 ], &q_1[ 0 ], &r_1[ 0 ], &p_2[ 0 ], &q_2[ 0 ], &r_2[ 0 ], &coplanar, &s[ 0 ], &t[ 0 ] ) )
		return eCutResult::Empty;

	if( coplanar == 1 )
	{
		// if(tri_tri_overlap_test_3d(&p_1[0], &q_1[0], &r_1[0], &p_2[0], &q_2[0], &r_2[0])) //TODO?
		return eCutResult::Coplanar;
	}

	if( s[ 0 ] == t[ 0 ] && s[ 1 ] == t[ 1 ] && s[ 2 ] == t[ 2 ] )
		return eCutResult::Empty;

	auto is_collinear = [ & ]( const Vector3& st, const Vector3& q12 )
	{
		Vector3 cross = st.cross( q12 );
		if( cross[ 0 ] <= SCALAR_ZERO && cross[ 1 ] <= SCALAR_ZERO && cross[ 2 ] <= SCALAR_ZERO )
			return true;
		return false;
	};

	Vector3 st( t[ 0 ] - s[ 0 ], t[ 1 ] - s[ 1 ], t[ 2 ] - s[ 2 ] );
	if( is_collinear( st, q1 - q2 ) )
		return eCutResult::Edge0;
	if( is_collinear( st, q2 - q3 ) )
		return eCutResult::Edge1;
	if( is_collinear( st, q3 - q1 ) )
		return eCutResult::Edge2;

	return eCutResult::Face;
}

floatTetWild::eCutResult floatTetWild::is_tri_tri_cutted_hint(
  const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& q1, const Vector3& q2, const Vector3& q3, eCutResult hint, const Logger* log )
{
	std::array<Scalar, 3> p_1 = { { 0, 0, 0 } }, q_1 = { { 1, 0, 0 } }, r_1 = { { 0, 1, 0 } };
	std::array<Scalar, 3> p_2 = { { 0, 0, 0 } }, q_2 = { { -1, -1, 0 } }, r_2 = { { 1, 1, 0 } };
	bool coplanar = false;
	std::array<Scalar, 3> s = { { 0, 0, 0 } }, t = { { 0, 0, 0 } };

	for( int j = 0; j < 3; j++ )
	{
		p_1[ j ] = p1[ j ];
		q_1[ j ] = p2[ j ];
		r_1[ j ] = p3[ j ];
		p_2[ j ] = q1[ j ];
		q_2[ j ] = q2[ j ];
		r_2[ j ] = q3[ j ];
	}

	if( hint == eCutResult::Coplanar )
	{
		int t = get_t( p1, p2, p3 );

		if( is_tri_tri_cutted_2d( { { to_2d( p1, t ), to_2d( p2, t ), to_2d( p3, t ) } }, { { to_2d( q1, t ), to_2d( q2, t ), to_2d( q3, t ) } } ) )
			return eCutResult::Coplanar;
		return eCutResult::Empty;
	}

	int result = triangleIntersectionTest( &p_1[ 0 ], &q_1[ 0 ], &r_1[ 0 ], &p_2[ 0 ], &q_2[ 0 ], &r_2[ 0 ], &coplanar, &s[ 0 ], &t[ 0 ] );
	if( nullptr != log )
		log->logDebug( ">>result = %i", result );

	if( result != 1 )
		return eCutResult::Empty;

	if( std::abs( s[ 0 ] - t[ 0 ] ) <= SCALAR_ZERO && std::abs( s[ 1 ] - t[ 1 ] ) <= SCALAR_ZERO && std::abs( s[ 2 ] - t[ 2 ] ) <= SCALAR_ZERO )
		return eCutResult::Empty;

	if( hint == eCutResult::Edge0 )
		return eCutResult::Edge0;
	if( hint == eCutResult::Edge1 )
		return eCutResult::Edge1;
	if( hint == eCutResult::Edge2 )
		return eCutResult::Edge2;

	return eCutResult::Face;
}

void floatTetWild::get_bbox_face( const Vector3& p0, const Vector3& p1, const Vector3& p2, __m256d& min, __m256d& max )
{
	__m256d i, ax;
	using AvxMath::loadDouble3;
	i = ax = loadDouble3( p0.data() );

	__m256d v = loadDouble3( p1.data() );
	i = _mm256_min_pd( i, v );
	ax = _mm256_max_pd( ax, v );

	v = loadDouble3( p2.data() );
	i = _mm256_min_pd( i, v );
	ax = _mm256_max_pd( ax, v );

	min = i;
	max = ax;
}

void floatTetWild::get_bbox_tet( const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3, __m256d& min, __m256d& max )
{
	__m256d i, ax;
	static_assert( sizeof( Vector3 ) == 3 * sizeof( double ) );
	using AvxMath::loadDouble3;

	i = ax = loadDouble3( p0.data() );

	__m256d v = loadDouble3( p1.data() );
	i = _mm256_min_pd( i, v );
	ax = _mm256_max_pd( ax, v );

	v = loadDouble3( p2.data() );
	i = _mm256_min_pd( i, v );
	ax = _mm256_max_pd( ax, v );

	v = loadDouble3( p3.data() );
	i = _mm256_min_pd( i, v );
	ax = _mm256_max_pd( ax, v );

	min = i;
	max = ax;
}

bool floatTetWild::is_bbox_intersected( const Vector3& min1, const Vector3& max1, const Vector3& min2, const Vector3& max2 )
{
	for( int j = 0; j < 3; j++ )
	{
		if( min1[ j ] > max2[ j ] )
			return false;
		if( max1[ j ] < min2[ j ] )
			return false;
	}
	return true;
}

bool floatTetWild::is_bbox_intersected( __m256d min1, __m256d max1, __m256d min2, __m256d max2 )
{
	// Compare for min1 > max2
	__m256d cmp1 = _mm256_cmp_pd( min1, max2, _CMP_GT_OQ );
	// Compare for max1 < min2
	__m256d cmp2 = _mm256_cmp_pd( max1, min2, _CMP_LT_OQ );
	// Combine with bitwise OR
	__m256d cmp = _mm256_or_pd( cmp1, cmp2 );
	// Move bitmap to a scalar register
	uint32_t mask = (uint32_t)_mm256_movemask_pd( cmp );
	// Produce the result
	return 0 == ( mask & 0b111 );
}

bool floatTetWild::is_tri_inside_tet( const std::array<Vector3, 3>& ps, const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t )
{
	int cnt_pos = 0;
	int cnt_neg = 0;

	for( int i = 0; i < 3; i++ )
	{
		eOrientation ori = Predicates::orient_3d( ps[ i ], p1t, p2t, p3t );
		if( ori == eOrientation::Positive )
			cnt_pos++;
		else if( ori == eOrientation::Negative )
			cnt_neg++;

		ori = Predicates::orient_3d( p0t, ps[ i ], p2t, p3t );
		if( ori == eOrientation::Positive )
			cnt_pos++;
		else if( ori == eOrientation::Negative )
			cnt_neg++;

		ori = Predicates::orient_3d( p0t, p1t, ps[ i ], p3t );
		if( ori == eOrientation::Positive )
			cnt_pos++;
		else if( ori == eOrientation::Negative )
			cnt_neg++;

		ori = Predicates::orient_3d( p0t, p1t, p2t, ps[ i ] );
		if( ori == eOrientation::Positive )
			cnt_pos++;
		else if( ori == eOrientation::Negative )
			cnt_neg++;
	}

	if( cnt_pos == 0 || cnt_neg == 0 )
		return true;

	return false;
}

bool floatTetWild::is_point_inside_tet( const Vector3& p, const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t )
{  /// inside or on
	int cnt_pos = 0;
	int cnt_neg = 0;

	eOrientation ori = Predicates::orient_3d( p, p1t, p2t, p3t );
	if( ori == eOrientation::Positive )
		cnt_pos++;
	else if( ori == eOrientation::Negative )
		cnt_neg++;

	ori = Predicates::orient_3d( p0t, p, p2t, p3t );
	if( ori == eOrientation::Positive )
		cnt_pos++;
	else if( ori == eOrientation::Negative )
		cnt_neg++;

	ori = Predicates::orient_3d( p0t, p1t, p, p3t );
	if( ori == eOrientation::Positive )
		cnt_pos++;
	else if( ori == eOrientation::Negative )
		cnt_neg++;

	ori = Predicates::orient_3d( p0t, p1t, p2t, p );
	if( ori == eOrientation::Positive )
		cnt_pos++;
	else if( ori == eOrientation::Negative )
		cnt_neg++;

	if( cnt_pos == 0 || cnt_neg == 0 )
		return true;

	return false;
}