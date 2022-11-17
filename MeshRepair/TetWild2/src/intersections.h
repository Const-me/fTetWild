// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
#pragma once
#include "Mesh.h"
#include "../external/Predicates.h"

namespace floatTetWild
{
	enum struct eCutResult : uint8_t
	{
		Edge0 = 0,
		Edge1 = 1,
		Edge2 = 2,
		Face = 3,
		Coplanar = 4,
		Empty = 0xFF
	};

	eCutResult is_tri_tri_cutted( const std::array<Vector3, 3>& f_tri, const std::array<Vector3, 3>& f_tet, const std::array<eOrientation, 3>& oris_tri );

	Scalar seg_seg_squared_dist_3d( const std::array<Vector3, 2>& s1, const std::array<Vector3, 2>& s2 );

	Scalar p_seg_squared_dist_3d( const Vector3& p, const Vector3& a, const Vector3& b );
	Scalar p_line_squared_dist_3d( const Vector3& p, const Vector3& a, const Vector3& b );

	bool is_p_inside_tri_2d( const Vector2& p, const std::array<Vector2, 3>& tri );
	bool is_seg_tri_cutted_2d( const std::array<Vector2, 2>& seg, const std::array<Vector2, 3>& tri );
	bool is_tri_tri_cutted_2d( const std::array<Vector2, 3>& p_tet, const std::array<Vector2, 3>& p_tri );

	bool seg_seg_intersection_2d( const std::array<Vector2, 2>& seg1, const std::array<Vector2, 2>& seg2, Scalar& t2 );
	bool seg_line_intersection_2d( const std::array<Vector2, 2>& seg, const std::array<Vector2, 2>& line, Scalar& t_seg );
	bool seg_plane_intersection( const Vector3& p1, const Vector3& p2, const Vector3& a, const Vector3& n, Vector3& p, Scalar& d1 );

	int get_t( const Vector3& p0, const Vector3& p1, const Vector3& p2 );
	int get_t( const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, int idx );

	Vector2 to_2d( const Vector3& p, int t );
	Vector2 to_2d( const Vector3& p, const Vector3& n, const Vector3& pp, int t );

	bool is_crossing( eOrientation s1, eOrientation s2 );

	eCutResult is_tri_tri_cutted( const Vector3& p1, const Vector3& p2, const Vector3& p3,	// cutting tri
	  const Vector3& q1, const Vector3& q2, const Vector3& q3 );						  // face of tet
	eCutResult is_tri_tri_cutted_hint( const Vector3& p1, const Vector3& p2, const Vector3& p3,	 // cutting tri
	  const Vector3& q1, const Vector3& q2, const Vector3& q3, eCutResult hint,
	  const Logger* log = nullptr );  // face of tet

	void get_bbox_face( const Vector3& p0, const Vector3& p1, const Vector3& p2, __m256d& min, __m256d& max );
	void get_bbox_tet( const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3, __m256d& min, __m256d& max );

	bool is_bbox_intersected( const Vector3& min1, const Vector3& max1, const Vector3& min2, const Vector3& max2 );
	bool is_bbox_intersected( __m256d min1, __m256d max1, __m256d min2, __m256d max2 );

	bool is_tri_inside_tet( const std::array<Vector3, 3>& ps, const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t );
	bool is_point_inside_tet( const Vector3& p, const Vector3& p0t, const Vector3& p1t, const Vector3& p2t, const Vector3& p3t );
}  // namespace floatTetWild