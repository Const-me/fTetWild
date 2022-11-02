#pragma once
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/predicates.h>
#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_partition.h>
#include <geogram/points/nn_search.h>
#include <geogram/mesh/mesh_AABB.h>

namespace GEO2
{
	using GEO::Box;
	using GEO::coord_index_t;
	using GEO::index_t;

	using GEO::geo_sgn;
	using GEO::geo_sqr;

	using GEO::NO_FACET;
	using GEO::vec3;
	using GEO::vector;
	namespace PCK
	{
		using GEO::PCK::orient_3d;
	}
	using GEO::bbox_diagonal;
	using GEO::Delaunay;
	using GEO::Delaunay_var;
	using GEO::MESH_ORDER_MORTON;

	using GEO::cross;
	using GEO::distance;
	using GEO::distance2;
	using GEO::dot;
	using GEO::length;
	using GEO::normalize;

	using GEO::Attribute;
	using GEO::mesh_partition;
	using GEO::MESH_PARTITION_HILBERT;
	using GEO::MeshCellsAABB;
	using GEO::NearestNeighborSearch;
	using GEO::NearestNeighborSearch_var;

	namespace Geom
	{
		using GEO::Geom::distance2;
		using GEO::Geom::point_triangle_squared_distance;
		using GEO::Geom::tetra_signed_volume;
	}  // namespace Geom

	using TetraMesh = GEO::Mesh;

	class Mesh
	{
		GEO::Mesh m;

	  public:
		class FaceCollection
		{
			GEO::MeshFacets& mesh;

		  public:
			FaceCollection( GEO::MeshFacets& orig )
				: mesh( orig )
			{
			}

			index_t nb() const { return mesh.nb(); }
			index_t nb_vertices( index_t f ) const { return mesh.nb_vertices( f ); }
			index_t corners_begin( index_t f ) const { return mesh.corners_begin( f ); }
			index_t corners_end( index_t f ) const { return mesh.corners_end( f ); }

			void clear( bool keep_attributes = true, bool keep_memory = false ) { mesh.clear( keep_attributes, keep_memory ); }
			index_t create_triangles( index_t nb ) { return mesh.create_triangles( nb ); }
			void set_vertex( index_t f, index_t lv, index_t v ) { mesh.set_vertex( f, lv, v ); }
			index_t vertex( index_t f, index_t lv ) const { return mesh.vertex( f, lv ); }

			void assign_triangle_mesh( vector<index_t>& triangles, bool steal_args ) { mesh.assign_triangle_mesh( triangles, steal_args ); }
		};
		FaceCollection facets;

		class FacetCornersCollection
		{
			GEO::MeshFacetCornersStore& mesh;

		  public:
			FacetCornersCollection( GEO::MeshFacetCornersStore& orig )
				: mesh( orig )
			{
			}

			index_t vertex( index_t c ) const { return mesh.vertex( c ); }
		};
		FacetCornersCollection facet_corners;

		class VertexCollection
		{
			GEO::MeshVertices& mesh;

		  public:
			VertexCollection( GEO::MeshVertices& orig )
				: mesh( orig )
			{
			}

			index_t nb() const { return mesh.nb(); }

			double* point_ptr( index_t v ) { return mesh.point_ptr( v ); }
			const double* point_ptr( index_t v ) const { return mesh.point_ptr( v ); }

			vec3& point( index_t v ) { return mesh.point( v ); }
			const vec3& point( index_t v ) const { return mesh.point( v ); }

			void clear( bool keep_attributes = true, bool keep_memory = false ) { mesh.clear( keep_attributes, keep_memory ); }
			index_t create_vertices( index_t nb ) { return mesh.create_vertices( nb ); }

			void assign_points( vector<double>& points, index_t dim, bool steal_arg ) 
			{
				if( 3 != dim )
					__debugbreak();
				mesh.assign_points( points, 3, steal_arg );
			}
		};
		VertexCollection vertices;

		Mesh()
			: facets( m.facets )
			, facet_corners( m.facet_corners )
			, vertices( m.vertices )
		{
		}

		const vec3& getVertex( index_t v ) const { return GEO::Geom::mesh_vertex( m, v ); }
		void clear( bool keep_attributes = true, bool keep_memory = false ) { m.clear( keep_attributes, keep_memory ); }
		double boxDiagonal() const { return bbox_diagonal( m ); }

		void reorderMorton() { mesh_reorder( m, GEO2::MESH_ORDER_MORTON ); }
	};

	namespace Geom
	{
		inline const vec3& mesh_vertex( const GEO2::Mesh& M, index_t v )
		{
			return M.getVertex( v );
		}
	}  // namespace Geom

	inline double bbox_diagonal( const GEO2::Mesh& M )
	{
		return M.boxDiagonal();
	}
}  // namespace GEO2