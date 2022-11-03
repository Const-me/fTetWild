#pragma once
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_reorder.h>
#include "../../src/Types.hpp"

// A wrapper around geogram's mesh
class TriMeshWrapper
{
	GEO::Mesh m;

  public:
	using index_t = GEO::index_t;
	using vec3 = GEO::vec3;

	class FaceCollection
	{
		GEO::MeshFacets& mesh;

	  public:
		FaceCollection( GEO::MeshFacets& orig )
			: mesh( orig )
		{
		}

		index_t nb() const
		{
			return mesh.nb();
		}
		index_t nb_vertices( index_t f ) const
		{
			return mesh.nb_vertices( f );
		}
		index_t corners_begin( index_t f ) const
		{
			return mesh.corners_begin( f );
		}
		index_t corners_end( index_t f ) const
		{
			return mesh.corners_end( f );
		}

		void clear( bool keep_attributes = true, bool keep_memory = false )
		{
			mesh.clear( keep_attributes, keep_memory );
		}
		index_t create_triangles( index_t nb )
		{
			return mesh.create_triangles( nb );
		}
		void set_vertex( index_t f, index_t lv, index_t v )
		{
			mesh.set_vertex( f, lv, v );
		}
		index_t vertex( index_t f, index_t lv ) const
		{
			return mesh.vertex( f, lv );
		}

		void assign_triangle_mesh( GEO::vector<index_t>& triangles, bool steal_args )
		{
			mesh.assign_triangle_mesh( triangles, steal_args );
		}
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

		index_t vertex( index_t c ) const
		{
			return mesh.vertex( c );
		}
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

		index_t nb() const
		{
			return mesh.nb();
		}

		double* point_ptr( index_t v )
		{
			return mesh.point_ptr( v );
		}
		const double* point_ptr( index_t v ) const
		{
			return mesh.point_ptr( v );
		}

		vec3& point( index_t v )
		{
			return mesh.point( v );
		}
		const vec3& point( index_t v ) const
		{
			return mesh.point( v );
		}

		void clear( bool keep_attributes = true, bool keep_memory = false )
		{
			mesh.clear( keep_attributes, keep_memory );
		}
		index_t create_vertices( index_t nb )
		{
			return mesh.create_vertices( nb );
		}

		void assign_points( GEO::vector<double>& points, index_t dim, bool steal_arg )
		{
			if( 3 != dim )
				__debugbreak();
			mesh.assign_points( points, 3, steal_arg );
		}
	};
	VertexCollection vertices;

	TriMeshWrapper()
		: facets( m.facets )
		, facet_corners( m.facet_corners )
		, vertices( m.vertices )
	{
	}

	const vec3& getVertex( index_t v ) const
	{
		return GEO::Geom::mesh_vertex( m, v );
	}
	void clear( bool keep_attributes = true, bool keep_memory = false )
	{
		m.clear( keep_attributes, keep_memory );
	}
	double boxDiagonal() const
	{
		return bbox_diagonal( m );
	}

	void reorderMorton()
	{
		GEO::mesh_reorder( m, GEO::MESH_ORDER_MORTON );
	}

	uint32_t countTriangles() const
	{
		return m.facets.nb();
	}

	void getTriangleVertices( uint32_t tri, const GEO::vec3** p1, const GEO::vec3** p2, const GEO::vec3** p3 ) const
	{
		index_t c = m.facets.corners_begin( tri );
		*p1 = &GEO::Geom::mesh_vertex( m, m.facet_corners.vertex( c ) );
		++c;
		*p2 = &GEO::Geom::mesh_vertex( m, m.facet_corners.vertex( c ) );
		++c;
		*p3 = &GEO::Geom::mesh_vertex( m, m.facet_corners.vertex( c ) );
	}

	// Set vertex buffer of the mesh, upcasting coordinates to FP64
	HRESULT assignVertices( size_t count, const float* vb );
	// Set index buffer of the mesh
	HRESULT assignTriangles( size_t count, const uint32_t* ib );
	HRESULT copyData( std::vector<floatTetWild::Vector3>& vb, std::vector<floatTetWild::Vector3i>& ib ) const;

	const GEO::vec3& getFirstTriangleVertex( uint32_t tri ) const
	{
		index_t v = m.facet_corners.vertex( m.facets.corners_begin( tri ) );
		return GEO::Geom::mesh_vertex( m, v );
	}

	template<class Lambda>
	void generateVertices( uint32_t count, Lambda lambda )
	{
		m.vertices.clear();
		m.vertices.create_vertices( count );
		for( uint32_t i = 0; i < count; i++ )
			m.vertices.point( i ) = lambda( i );
	}

	template<class Lambda>
	void generateTriangles( uint32_t count, Lambda lambda )
	{
		m.facets.clear();
		m.facets.create_triangles( count );
		for( uint32_t i = 0; i < count; i++ ) 
		{
			const GEO::vec3i tri = lambda( i );
			m.facets.set_vertex( i, 0, tri.x );
			m.facets.set_vertex( i, 1, tri.y );
			m.facets.set_vertex( i, 2, tri.z );
		}
	}

	void clearMesh()
	{
		m.clear( false, true );
	}
};