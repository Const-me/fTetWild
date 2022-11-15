#include "stdafx.h"
#include "LocalOperations.h"
#include "AABBWrapper.h"
#include "TriangleInsertion.h"

void setSingleTriangle( GEO2::Mesh& m )
{
	m.generateVertices( 1, []( uint32_t ) { return GEO2::vec3( 0, 0, 0 ); } );
	m.generateTriangles( 1, []( uint32_t ) { return GEO2::vec3i { 0, 0, 0 }; } );
}

void generateEdgeTriangles( GEO2::Mesh& m, uint32_t count )
{
	m.generateTriangles( count,
	  []( uint32_t i )
	  {
		  GEO2::vec3i res;
		  res.x = res.y = (int)i * 2;
		  res.z = res.x + 1;
		  return res;
	  } );
}

inline GEO2::vec3 castVertex( const floatTetWild::Vector3& eigen )
{
	return GEO2::vec3 { eigen[ 0 ], eigen[ 1 ], eigen[ 2 ] };
}

void floatTetWild::AABBWrapper::init_b_mesh_and_tree(
  const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, Mesh& mesh )
{
	b_mesh.clearMesh();
	std::vector<std::vector<int>> conn_tris( input_vertices.size() );
	std::vector<std::array<int, 2>> all_edges;
	all_edges.reserve( input_faces.size() * 3 );
	for( int i = 0; i < input_faces.size(); i++ )
	{
		for( int j = 0; j < 3; j++ )
		{
			conn_tris[ input_faces[ i ][ j ] ].push_back( i );
			if( input_faces[ i ][ j ] < input_faces[ i ][ ( j + 1 ) % 3 ] )
				all_edges.push_back( { { input_faces[ i ][ j ], input_faces[ i ][ ( j + 1 ) % 3 ] } } );
			else
				all_edges.push_back( { { input_faces[ i ][ ( j + 1 ) % 3 ], input_faces[ i ][ j ] } } );
		}
	}
	vector_unique( all_edges );

	std::vector<std::pair<std::array<int, 2>, std::vector<int>>> _;
	std::vector<std::array<int, 2>> b_edges;
	std::vector<bool> _1;
	find_boundary_edges(
	  input_vertices, input_faces, std::vector<bool>( input_faces.size(), true ), std::vector<bool>( input_faces.size(), true ), _, _1, b_edges, mesh.logger() );

	if( b_edges.empty() )
		setSingleTriangle( b_mesh );
	else
	{
		b_mesh.generateVertices( b_edges.size() * 2,
		  [ & ]( uint32_t i )
		  {
			  const auto& edge = b_edges[ i / 2 ];
			  return castVertex( input_vertices[ edge[ i % 2 ] ] );
		  } );

		generateEdgeTriangles( b_mesh, b_edges.size() );
	}

	b_mesh.reorderMorton();
	b_tree = std::make_shared<MeshFacetsAABBWithEps>( b_mesh, mesh.facetRecursionStacks );

	if( b_edges.empty() )
		mesh.is_closed = true;
}

void floatTetWild::AABBWrapper::init_tmp_b_mesh_and_tree( const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces,
  const std::vector<std::array<int, 2>>& b_edges1, const Mesh& mesh, const std::vector<std::array<int, 2>>& b_edges2 )
{
	if( b_edges1.empty() && b_edges2.empty() )
		setSingleTriangle( tmp_b_mesh );
	else
	{
		tmp_b_mesh.generateVertices( ( b_edges1.size() + b_edges2.size() ) * 2,
		  [ & ]( uint32_t i )
		  {
			  const uint32_t edge = i / 2;
			  if( edge < b_edges1.size() )
			  {
				  const auto& e = b_edges1[ edge ];
				  return castVertex( input_vertices[ e[ i % 2 ] ] );
			  }
			  else
			  {
				  const auto& e = b_edges2[ edge - b_edges1.size() ];
				  return castVertex( mesh.tet_vertices[ e[ i % 2 ] ].pos );
			  }
		  } );

		generateEdgeTriangles( tmp_b_mesh, b_edges1.size() + b_edges2.size() );
	}
	tmp_b_mesh.reorderMorton();
	tmp_b_tree = std::make_shared<MeshFacetsAABBWithEps>( tmp_b_mesh, const_cast<FacetRecursionStacks&>( mesh.facetRecursionStacks ) );
}