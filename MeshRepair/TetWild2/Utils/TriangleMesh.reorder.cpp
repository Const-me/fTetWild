#include "stdafx.h"
#include "TriangleMesh.h"
#include <numeric>
#include "reorderUtils.h"
#include <geogram/basic/process.h>

// The implementation was copy-pasted from there:
// https://github.com/BrunoLevy/geogram/blob/main/src/lib/geogram/mesh/mesh_reorder.cpp
// BSD license allows that
namespace
{
	using GEO::vec3;
	using GEO::vec3i;
	using GEO::vector;

	/**
	 * \brief Splits a sequence into two ordered halves.
	 * \details The algorithm shuffles the sequence and
	 *  partitions its into two halves with the same number of elements
	 *  and such that the elements of the first half are smaller
	 *  than the elements of the second half.
	 * \param[in] begin an iterator to the first element
	 * \param[in] end an iterator one position past the last element
	 * \param[in] cmp the comparator object
	 * \return an iterator to the middle of the sequence that separates
	 *  the two halves
	 */
	template<class IT, class CMP>
	inline IT reorder_split( IT begin, IT end, CMP cmp )
	{
		if( begin >= end )
		{
			return begin;
		}
		IT middle = begin + ( end - begin ) / 2;
		std::nth_element( begin, middle, end, cmp );
		return middle;
	}

	template<int COORD, bool UP>
	struct Hilbert_vcmp
	{
		Hilbert_vcmp( const vector<vec3>& vertexBuffer )
			: vb( vertexBuffer )
		{
			static_assert( COORD >= 0 && COORD < 3, "Wrong template argument" );
		}

		bool operator()( uint32_t i1, uint32_t i2 )
		{
			double a = load( i1 );
			double b = load( i2 );
			if( UP )  // TODO: upgrade to C++/17, replace with if constexpr
				return a < b;
			else
				return a > b;
		}

	  private:
		const vector<vec3>& vb;

		double load( uint32_t idx ) const
		{
			const double* const ptr = (const double*)&vb[ idx ];
			return ptr[ COORD ];
		}
	};

	template<int COORD, bool UP>
	using Morton_vcmp = Hilbert_vcmp<COORD, true>;

	template<template<int COORD, bool UP> class CMP, class MESH>
	class HilbertSort3d
	{
		template<int COORDX, bool UPX, bool UPY, bool UPZ, class IT>
		static void sort( const MESH& M, IT begin, IT end, ptrdiff_t limit = 1 )
		{
			const int COORDY = ( COORDX + 1 ) % 3, COORDZ = ( COORDY + 1 ) % 3;
			if( end - begin <= limit )
				return;

			IT m0 = begin, m8 = end;
			IT m4 = reorder_split( m0, m8, CMP<COORDX, UPX>( M ) );
			IT m2 = reorder_split( m0, m4, CMP<COORDY, UPY>( M ) );
			IT m1 = reorder_split( m0, m2, CMP<COORDZ, UPZ>( M ) );
			IT m3 = reorder_split( m2, m4, CMP<COORDZ, !UPZ>( M ) );
			IT m6 = reorder_split( m4, m8, CMP<COORDY, !UPY>( M ) );
			IT m5 = reorder_split( m4, m6, CMP<COORDZ, UPZ>( M ) );
			IT m7 = reorder_split( m6, m8, CMP<COORDZ, !UPZ>( M ) );
			sort<COORDZ, UPZ, UPX, UPY>( M, m0, m1 );
			sort<COORDY, UPY, UPZ, UPX>( M, m1, m2 );
			sort<COORDY, UPY, UPZ, UPX>( M, m2, m3 );
			sort<COORDX, UPX, !UPY, !UPZ>( M, m3, m4 );
			sort<COORDX, UPX, !UPY, !UPZ>( M, m4, m5 );
			sort<COORDY, !UPY, UPZ, !UPX>( M, m5, m6 );
			sort<COORDY, !UPY, UPZ, !UPX>( M, m6, m7 );
			sort<COORDZ, !UPZ, !UPX, UPY>( M, m7, m8 );
		}

	  public:
		HilbertSort3d( const MESH& M, vector<uint32_t>::iterator b, vector<uint32_t>::iterator e, ptrdiff_t limit = 1 )
			: M_( M )
		{
			assert( e >= b );

			// If the sequence is smaller than the limit, skip it
			if( ( e - b ) <= limit )
				return;

			// If the sequence is smaller than 1024, use sequential sorting
			if( ( e - b ) < 1024 )
			{
				sort<0, false, false, false>( M_, b, e );
				return;
			}
#if 0
			// Parallel sorting (2 then 4 then 8 sorts in parallel)
			// Unfortunately we cannot access consts for template arguments in lambdas in all
			// compilers (gcc is OK but not MSVC) so I'm using (ugly) macros here...

#define COORDX 0
#define COORDY 1
#define COORDZ 2
#define UPX false
#define UPY false
#define UPZ false

			m0_ = b;
			m8_ = e;
			m4_ = reorder_split( m0_, m8_, CMP<COORDX, UPX>( M ) );

			parallel_for( [ this ]() { m2_ = reorder_split( m0_, m4_, CMP<COORDY, UPY>( M_ ) ); },
			  [ this ]() { m6_ = reorder_split( m4_, m8_, CMP<COORDY, !UPY>( M_ ) ); } );

			parallel_for( [ this ]() { m1_ = reorder_split( m0_, m2_, CMP<COORDZ, UPZ>( M_ ) ); },
			  [ this ]() { m3_ = reorder_split( m2_, m4_, CMP<COORDZ, !UPZ>( M_ ) ); },
			  [ this ]() { m5_ = reorder_split( m4_, m6_, CMP<COORDZ, UPZ>( M_ ) ); },
			  [ this ]() { m7_ = reorder_split( m6_, m8_, CMP<COORDZ, !UPZ>( M_ ) ); } );

			parallel_for( [ this ]() { sort<COORDZ, UPZ, UPX, UPY>( M_, m0_, m1_ ); }, [ this ]() { sort<COORDY, UPY, UPZ, UPX>( M_, m1_, m2_ ); },
			  [ this ]() { sort<COORDY, UPY, UPZ, UPX>( M_, m2_, m3_ ); }, [ this ]() { sort<COORDX, UPX, !UPY, !UPZ>( M_, m3_, m4_ ); },
			  [ this ]() { sort<COORDX, UPX, !UPY, !UPZ>( M_, m4_, m5_ ); }, [ this ]() { sort<COORDY, !UPY, UPZ, !UPX>( M_, m5_, m6_ ); },
			  [ this ]() { sort<COORDY, !UPY, UPZ, !UPX>( M_, m6_, m7_ ); }, [ this ]() { sort<COORDZ, !UPZ, !UPX, UPY>( M_, m7_, m8_ ); } );

#undef COORDX
#undef COORDY
#undef COORDZ
#undef UPX
#undef UPY
#undef UPZ
#else
			sort<0, false, false, false>( M_, b, e );
#endif
		}

	  private:
		const MESH& M_;
		vector<uint32_t>::iterator m0_, m1_, m2_, m3_, m4_, m5_, m6_, m7_, m8_;
	};

	void morton_vsort_3d( const vector<vec3>& vec, std::vector<uint32_t>& reordered )
	{
		reordered.resize( vec.size() );
		std::iota( reordered.begin(), reordered.end(), (uint32_t)0 );
		// HilbertSort3d<Hilbert_vcmp, vector<vec3>> sorter( vec, reordered.begin(), reordered.end() );
		HilbertSort3d<Morton_vcmp, vector<vec3>> sorter( vec, reordered.begin(), reordered.end() );
	}
}  // namespace

template<int COORD>
inline double TriangleMesh::triangleCenterX3( uint32_t idxTri ) const
{
	static_assert( COORD >= 0 && COORD < 3, "Invalid template parameter" );
	const vec3i& tri = triangles[ idxTri ];

	uint32_t v0 = (uint32_t)tri.x;
	uint32_t v1 = (uint32_t)tri.y;
	uint32_t v2 = (uint32_t)tri.z;
	assert( v0 < vertices.size() && v1 < vertices.size() && v2 < vertices.size() );

	const double* p0 = (const double*)&vertices[ v0 ];
	const double* p1 = (const double*)&vertices[ v1 ];
	const double* p2 = (const double*)&vertices[ v2 ];
	return p0[ COORD ] + p1[ COORD ] + p2[ COORD ];
}

namespace
{
	template<int COORD, bool UP>
	struct Hilbert_fcmp
	{
		Hilbert_fcmp( const TriangleMesh& m )
			: mesh( m )
		{
			static_assert( COORD >= 0 && COORD < 3, "Wrong template argument" );
		}

		bool operator()( uint32_t i1, uint32_t i2 )
		{
			double a = load( i1 );
			double b = load( i2 );
			if( UP )  // TODO: upgrade to C++/17, replace with if constexpr
				return a < b;
			else
				return a > b;
		}

	  private:
		const TriangleMesh& mesh;

		double load( uint32_t idx ) const
		{
			return mesh.triangleCenterX3<COORD>( idx );
		}
	};

	template<int COORD, bool UP>
	using Morton_fcmp = Hilbert_fcmp<COORD, true>;

	void morton_fsort_3d( const TriangleMesh& mesh, std::vector<uint32_t>& reordered )
	{
		reordered.resize( mesh.countTriangles() );
		std::iota( reordered.begin(), reordered.end(), (uint32_t)0 );

		HilbertSort3d<Morton_fcmp, TriangleMesh> sorter( mesh, reordered.begin(), reordered.end() );
	}
}  // namespace

void TriangleMesh::reorderMorton()
{
	std::vector<uint32_t> order;

	// Step #1, reorder vertices
	morton_vsort_3d( vertices, order );
	reorderVector( vertices, order );
	remapIndices( triangles, order );

	// Step #2, reorder triangles
	morton_fsort_3d( *this, order );
	reorderVector( triangles, order );
}