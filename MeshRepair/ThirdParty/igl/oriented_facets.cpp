// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "oriented_facets.h"

template<typename DerivedF, typename DerivedE>
IGL_INLINE void igl::oriented_facets( const Eigen::MatrixBase<DerivedF>& F, Eigen::PlainObjectBase<DerivedE>& E )
{
	const auto rows = F.rows();
	E.resize( rows * F.cols(), F.cols() - 1 );
	using EScalar = typename DerivedE::Scalar;
	switch( F.cols() )
	{
	case 3:
		E.block( 0 * rows, 0, rows, 1 ) = F.col( 1 ).template cast<EScalar>();
		E.block( 0 * rows, 1, rows, 1 ) = F.col( 2 ).template cast<EScalar>();
		E.block( 1 * rows, 0, rows, 1 ) = F.col( 2 ).template cast<EScalar>();
		E.block( 1 * rows, 1, rows, 1 ) = F.col( 0 ).template cast<EScalar>();
		E.block( 2 * rows, 0, rows, 1 ) = F.col( 0 ).template cast<EScalar>();
		E.block( 2 * rows, 1, rows, 1 ) = F.col( 1 ).template cast<EScalar>();
		return;
	default:
		throw std::exception( "igl::oriented_facets, unexpected matrix size" );
	}
}