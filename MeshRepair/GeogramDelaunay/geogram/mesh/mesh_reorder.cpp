#include "mesh_reorder.h"
#include <random>
#include <assert.h>

namespace GEO
{
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

	/**
	 * \brief Used by VertexMesh.
	 * \details Exposes an interface compatible with the requirement
	 * of Hilbert sort templates for a raw array of vertices.
	 */
	class VertexArray
	{
	  public:
		/**
		 * \brief Constructs a new VertexArray.
		 * \param[in] base address of the points
		 * \param[in] stride number of doubles between
		 *  two consecutive points
		 */
		VertexArray( index_t nb_vertices, const double* base, index_t stride )
			: base_( base )
			, stride_( stride )
		{
			nb_vertices_ = nb_vertices;
		}

		/**
		 * \brief Gets a vertex by its index.
		 * \param[in] i the index of the point
		 * \return a const pointer to the coordinates of the vertex
		 */
		const double* point_ptr( index_t i ) const
		{
			geo_debug_assert( i < nb_vertices_ );
			return base_ + i * stride_;
		}

	  private:
		const double* base_;
		index_t stride_;
		index_t nb_vertices_;
	};

	class VertexMesh
	{
	  public:
		/**
		 * \brief Constructs a new VertexMesh.
		 * \param[in] base address of the points
		 * \param[in] stride number of doubles between
		 *  two consecutive points
		 */
		VertexMesh( index_t nb_vertices, const double* base, index_t stride )
			: vertices( nb_vertices, base, stride )
		{
		}
		VertexArray vertices;
	};

	constexpr int COORDX = 0, COORDY = 1, COORDZ = 2;
	constexpr bool UPX = false, UPY = false, UPZ = false;

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

		void threadProc( int i )
		{
			switch( i )
			{
			case 0:
				m2_ = reorder_split( m0_, m4_, CMP<COORDY, UPY>( M_ ) );
				break;
			case 1:
				m6_ = reorder_split( m4_, m8_, CMP<COORDY, !UPY>( M_ ) );
				break;
			case 10:
				m1_ = reorder_split( m0_, m2_, CMP<COORDZ, UPZ>( M_ ) );
				break;
			case 11:
				m3_ = reorder_split( m2_, m4_, CMP<COORDZ, !UPZ>( M_ ) );
				break;
			case 12:
				m5_ = reorder_split( m4_, m6_, CMP<COORDZ, UPZ>( M_ ) );
				break;
			case 13:
				m7_ = reorder_split( m6_, m8_, CMP<COORDZ, !UPZ>( M_ ) );
				break;
			case 20:
				sort<COORDZ, UPZ, UPX, UPY>( M_, m0_, m1_ );
				break;
			case 21:
				sort<COORDY, UPY, UPZ, UPX>( M_, m1_, m2_ );
				break;
			case 22:
				sort<COORDY, UPY, UPZ, UPX>( M_, m2_, m3_ );
				break;
			case 23:
				sort<COORDX, UPX, !UPY, !UPZ>( M_, m3_, m4_ );
				break;
			case 24:
				sort<COORDX, UPX, !UPY, !UPZ>( M_, m4_, m5_ );
				break;
			case 25:
				sort<COORDY, !UPY, UPZ, !UPX>( M_, m5_, m6_ );
				break;
			case 26:
				sort<COORDY, !UPY, UPZ, !UPX>( M_, m6_, m7_ );
				break;
			case 27:
				sort<COORDZ, !UPZ, !UPX, UPY>( M_, m7_, m8_ );
				break;
			default:
				__debugbreak();
				break;
			}
		}

		void sortInParallel()
		{
			// computes m2, m6 in parallel
#pragma omp parallel for
			for( int i = 0; i < 2; i++ )
				threadProc( i );

				// computes m1,m3,m5,m7 in parallel
#pragma omp parallel for
			for( int i = 10; i < 14; i++ )
				threadProc( i );

				// sorts the 8 subsets in parallel
#pragma omp parallel for
			for( int i = 20; i < 28; i++ )
				threadProc( i );
		}

	  public:
		HilbertSort3d( const MESH& M, std::vector<uint32_t>::iterator b, std::vector<uint32_t>::iterator e, ptrdiff_t limit = 1 )
			: M_( M )
		{
			assert( e >= b );

			// If the sequence is smaller than the limit, skip it
			if( ( e - b ) <= limit )
				return;

			// If the sequence is smaller than 1024, use sequential sorting
			if( ( e - b ) < 1024 )
				sort<0, false, false, false>( M_, b, e );
			else
			{
				m0_ = b;
				m8_ = e;
				m4_ = reorder_split( m0_, m8_, CMP<COORDX, UPX>( M ) );
				sortInParallel();
			}
		}

	  private:
		const MESH& M_;
		std::vector<uint32_t>::iterator m0_, m1_, m2_, m3_, m4_, m5_, m6_, m7_, m8_;
	};

	template<int COORD, bool UP>
	struct Hilbert_vcmp
	{
		Hilbert_vcmp( const VertexMesh& vertexBuffer )
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
		const VertexMesh& vb;

		double load( uint32_t idx ) const
		{
			const double* const ptr = vb.vertices.point_ptr( idx );
			return ptr[ COORD ];
		}
	};

	void compute_BRIO_order_recursive( index_t nb_vertices, const double* vertices, index_t dimension, index_t stride, vector<index_t>& sorted_indices,
	  vector<index_t>::iterator b, vector<index_t>::iterator e, index_t threshold, double ratio, index_t& depth, vector<index_t>* levels )
	{
		geo_debug_assert( e > b );

		vector<index_t>::iterator m = b;
		if( index_t( e - b ) > threshold )
		{
			++depth;
			m = b + int( double( e - b ) * ratio );
			compute_BRIO_order_recursive( nb_vertices, vertices, dimension, stride, sorted_indices, b, m, threshold, ratio, depth, levels );
		}

		VertexMesh M( nb_vertices, vertices, stride );
		if( dimension == 3 )
		{
			HilbertSort3d<Hilbert_vcmp, VertexMesh>( M, m, e );
		}
		/* else if( dimension == 2 )
		{
			HilbertSort2d<Hilbert_vcmp, VertexMesh>( M, m, e );
		} */
		else
		{
			geo_assert_not_reached;
		}

		if( levels != nullptr )
			levels->push_back( index_t( e - sorted_indices.begin() ) );
	}

	void compute_BRIO_order( index_t nb_vertices, const double* vertices, vector<index_t>& sorted_indices, index_t dimension, index_t stride, index_t threshold,
	  double ratio, vector<index_t>* levels )
	{
		if( levels != nullptr )
		{
			levels->clear();
			levels->push_back( 0 );
		}
		index_t depth = 0;
		sorted_indices.resize( nb_vertices );
		for( index_t i = 0; i < nb_vertices; ++i )
			sorted_indices[ i ] = i;

		// The next three lines replace the following commented-out line
		//(random_shuffle is deprecated in C++17, and they call this
		//  progess...)
		// std::random_shuffle(sorted_indices.begin(), sorted_indices.end());
		std::random_device rng;
		std::mt19937 urng( rng() );
		std::shuffle( sorted_indices.begin(), sorted_indices.end(), urng );

		compute_BRIO_order_recursive(
		  nb_vertices, vertices, dimension, stride, sorted_indices, sorted_indices.begin(), sorted_indices.end(), threshold, ratio, depth, levels );
	}
}  // namespace GEO