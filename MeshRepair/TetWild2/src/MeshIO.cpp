// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Teseo Schneider <teseo.schneider@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include "stdafx.h"
#include "MeshIO.hpp"
#include "MeshImprovement.h"
#include "../external/MshSaver.h"
#include "Logger.hpp"

#include <igl/Timer.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <igl/write_triangle_mesh.h>

#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/numerics/predicates.h>

#include <numeric>

namespace floatTetWild
{
	namespace
	{

		void extract_volume_mesh(
		  const Mesh& mesh, const std::function<bool( int )>& skip_tet, const std::function<bool( int )>& skip_vertex, MatrixXs& V, Eigen::MatrixXi& T )
		{
			const auto& points = mesh.tet_vertices;
			const auto& tets = mesh.tets;

			V.resize( points.size(), 3 );
			T.resize( tets.size(), 4 );

			size_t index = 0;
			std::vector<int> old_2_new( points.size(), -1 );
			for( size_t i = 0; i < points.size(); ++i )
			{
				if( skip_vertex( i ) )
				{
					continue;
				}
				old_2_new[ i ] = index;
				V.row( index ) = points[ i ].pos.transpose();
				++index;
			}

			V.conservativeResize( index, 3 );

			index = 0;
			for( size_t i = 0; i < tets.size(); ++i )
			{
				if( skip_tet( i ) )
					continue;
				for( int j = 0; j < 4; j++ )
				{
					T( index, j ) = old_2_new[ tets[ i ][ j ] ];
				}
				++index;
			}
			T.conservativeResize( index, 4 );
		}

		void extract_surface_mesh( const Mesh& mesh, const std::function<bool( int )>& skip_tet, const std::function<bool( int )>& skip_vertex,
		  Eigen::Matrix<Scalar, Eigen::Dynamic, 3>& VS, Eigen::Matrix<int, Eigen::Dynamic, 3>& FS )
		{
			MatrixXs VT;
			Eigen::MatrixXi TT;
			extract_volume_mesh( mesh, skip_tet, skip_vertex, VT, TT );

			Eigen::VectorXi I;
			igl::boundary_facets( TT, FS );
			igl::remove_unreferenced( VT, FS, VS, FS, I );
		}

		void write_mesh_aux( const std::string& path, const Mesh& mesh, const std::vector<int>& t_ids, const std::vector<Scalar>& color, const bool binary,
		  const bool separate_components, const std::function<bool( int )>& skip_tet, const std::function<bool( int )>& skip_vertex )
		{
			std::string output_format = path.substr( path.size() - 4, 4 );

			if( output_format == "mesh" )
			{
				if( binary )
					logger().warn( "Only non binary mesh supported, ignoring" );

				std::ofstream f( path );
				f.precision( std::numeric_limits<Scalar>::digits10 + 1 );

				f << "MeshVersionFormatted 1" << std::endl;
				f << "Dimension 3" << std::endl;

				int cnt_v = 0;
				std::map<int, int> old_2_new;
				for( int i = 0; i < mesh.tet_vertices.size(); i++ )
				{
					if( !skip_vertex( i ) )
					{
						old_2_new[ i ] = cnt_v;
						cnt_v++;
					}
				}
				int cnt_t = 0;
				for( const int i : t_ids )
				{
					if( !skip_tet( i ) )
						cnt_t++;
				}

				f << "Vertices" << std::endl << cnt_v << std::endl;

				for( size_t i = 0; i < mesh.tet_vertices.size(); i++ )
				{
					if( skip_vertex( i ) )
						continue;
					f << mesh.tet_vertices[ i ][ 0 ] << " " << mesh.tet_vertices[ i ][ 1 ] << " " << mesh.tet_vertices[ i ][ 2 ] << " " << 0 << std::endl;
				}

				f << "Triangles" << std::endl << 0 << std::endl;
				f << "Tetrahedra" << std::endl << cnt_t << std::endl;
				const std::array<int, 4> new_indices = { { 0, 1, 3, 2 } };

				for( const int i : t_ids )
				{
					if( skip_tet( i ) )
						continue;
					for( int j = 0; j < 4; j++ )
					{
						f << old_2_new[ mesh.tets[ i ][ new_indices[ j ] ] ] + 1 << " ";
					}
					f << 0 << std::endl;
				}

				f << "End";
				f.close();
			}
			else
			{
				assert( color.empty() || color.size() == mesh.tet_vertices.size() || color.size() == mesh.tets.size() );

				PyMesh::MshSaver mesh_saver( path, binary );

				std::map<int, int> old_2_new;
				int cnt_v = 0;
				for( int i = 0; i < mesh.tet_vertices.size(); i++ )
				{
					if( !skip_vertex( i ) )
					{
						old_2_new[ i ] = cnt_v;
						cnt_v++;
					}
				}
				int cnt_t = 0;
				for( const int i : t_ids )
				{
					if( !skip_tet( i ) )
						cnt_t++;
				}
				PyMesh::VectorF V_flat( cnt_v * 3 );
				PyMesh::VectorI T_flat( cnt_t * 4 );
				PyMesh::VectorI C_flat;

				if( separate_components )
					C_flat.resize( cnt_t );

				size_t index = 0;
				for( size_t i = 0; i < mesh.tet_vertices.size(); ++i )
				{
					if( skip_vertex( i ) )
						continue;
					for( int j = 0; j < 3; j++ )
						V_flat[ index * 3 + j ] = mesh.tet_vertices[ i ][ j ];
					index++;
				}

				index = 0;
				for( const int i : t_ids )
				{
					if( skip_tet( i ) )
						continue;
					T_flat[ index * 4 + 0 ] = old_2_new[ mesh.tets[ i ][ 0 ] ];
					T_flat[ index * 4 + 1 ] = old_2_new[ mesh.tets[ i ][ 1 ] ];
					T_flat[ index * 4 + 2 ] = old_2_new[ mesh.tets[ i ][ 3 ] ];
					T_flat[ index * 4 + 3 ] = old_2_new[ mesh.tets[ i ][ 2 ] ];

					if( separate_components )
						C_flat[ index ] = mesh.tets[ i ].scalar;

					index++;
				}

				mesh_saver.save_mesh( V_flat, T_flat, C_flat, 3, mesh_saver.TET );

				if( color.size() == mesh.tets.size() )
				{
					PyMesh::VectorF color_flat( cnt_t );
					index = 0;
					for( const int i : t_ids )
					{
						if( skip_tet( i ) )
							continue;
						color_flat[ index++ ] = color[ i ];
					}
					mesh_saver.save_elem_scalar_field( "color", color_flat );
				}
				else if( color.size() == mesh.tet_vertices.size() )
				{
					PyMesh::VectorF color_flat( cnt_v );
					index = 0;
					for( int i = 0; i < mesh.tet_vertices.size(); i++ )
					{
						if( skip_vertex( i ) )
							continue;
						color_flat[ index++ ] = color[ i ];
					}
					mesh_saver.save_scalar_field( "color", color_flat );
				}
			}
		}
	}  // namespace

	bool MeshIO::load_mesh( const std::string& path, std::vector<Vector3>& points, std::vector<Vector3i>& faces, GEO2::Mesh& input, std::vector<int>& flags )
	{
		__debugbreak();
		return false;
	}

	bool MeshIO::load_mesh( const std::string& path, std::vector<Vector3>& points, std::vector<Vector3i>& faces, GEO2::Mesh& input, std::vector<int>& flags,
	  std::vector<double>& epsr_flags )
	{
		__debugbreak();
		return false;
	}

	void MeshIO::load_mesh( std::vector<Vector3>& points, std::vector<Vector3i>& faces, GEO2::Mesh& input, std::vector<int>& flags )
	{
		__debugbreak();
	}

	void MeshIO::load_mesh(
	  std::vector<Vector3>& points, std::vector<Vector3i>& faces, GEO2::Mesh& input, std::vector<int>& flags, std::vector<double>& epsr_flags )
	{
		__debugbreak();
	}

	void MeshIO::write_mesh(
	  const std::string& path, const Mesh& mesh, const std::vector<int>& t_ids, const bool only_interior, const bool binary, const bool separate_components )
	{
		logger().info( "Writing mesh to {}...", path );
		igl::Timer timer;
		timer.start();

		if( only_interior )
		{
			const auto skip_tet = [ &mesh ]( const int i ) { return mesh.tets[ i ].is_outside; };
			const auto skip_vertex = [ &mesh ]( const int i ) { return mesh.tet_vertices[ i ].is_outside; };
			write_mesh_aux( path, mesh, t_ids, std::vector<Scalar>(), binary, separate_components, skip_tet, skip_vertex );
		}
		else
		{
			timer.start();
			const auto skip_tet = [ &mesh ]( const int i ) { return mesh.tets[ i ].is_removed; };
			const auto skip_vertex = [ &mesh ]( const int i ) { return mesh.tet_vertices[ i ].is_removed; };
			write_mesh_aux( path, mesh, t_ids, std::vector<Scalar>(), binary, separate_components, skip_tet, skip_vertex );
		}

		timer.stop();
		logger().info( " took {}s", timer.getElapsedTime() );
	}

	void MeshIO::write_mesh(
	  const std::string& path, const Mesh& mesh, const bool only_interior, const std::vector<Scalar>& color, const bool binary, const bool separate_components )
	{
		logger().info( "Writing mesh to {}...", path );
		igl::Timer timer;
		timer.start();

		std::vector<int> t_ids( mesh.tets.size() );
		std::iota( std::begin( t_ids ), std::end( t_ids ), 0 );	 // Fill with 0, 1, ..., n.

		if( only_interior )
		{
			const auto skip_tet = [ &mesh ]( const int i ) { return mesh.tets[ i ].is_outside; };
			const auto skip_vertex = [ &mesh ]( const int i ) { return mesh.tet_vertices[ i ].is_outside; };
			write_mesh_aux( path, mesh, t_ids, color, binary, separate_components, skip_tet, skip_vertex );
		}
		else
		{
			const auto skip_tet = [ &mesh ]( const int i ) { return mesh.tets[ i ].is_removed; };
			const auto skip_vertex = [ &mesh ]( const int i ) { return mesh.tet_vertices[ i ].is_removed; };
			write_mesh_aux( path, mesh, t_ids, color, binary, separate_components, skip_tet, skip_vertex );
		}

		timer.stop();
		logger().info( " took {}s", timer.getElapsedTime() );
	}

	void MeshIO::write_surface_mesh( const std::string& path, const Mesh& mesh, const bool only_interior )
	{
		logger().debug( "Extracting and writing surface to {}...", path );
		igl::Timer timer;
		timer.start();

		Eigen::Matrix<Scalar, Eigen::Dynamic, 3> V_sf;
		Eigen::Matrix<int, Eigen::Dynamic, 3> F_sf;
		if( only_interior )
		{
			const auto skip_tet = [ &mesh ]( const int i ) { return mesh.tets[ i ].is_outside; };
			const auto skip_vertex = [ &mesh ]( const int i ) { return mesh.tet_vertices[ i ].is_outside; };
			extract_surface_mesh( mesh, skip_tet, skip_vertex, V_sf, F_sf );
		}
		else
		{
			const auto skip_tet = [ &mesh ]( const int i ) { return mesh.tets[ i ].is_removed; };
			const auto skip_vertex = [ &mesh ]( const int i ) { return mesh.tet_vertices[ i ].is_removed; };
			extract_surface_mesh( mesh, skip_tet, skip_vertex, V_sf, F_sf );
		}

		igl::write_triangle_mesh( path, V_sf, F_sf );

		timer.stop();
		logger().info( " took {}s", timer.getElapsedTime() );
	}

	void MeshIO::extract_volume_mesh( const Mesh& mesh, MatrixXs& V, Eigen::MatrixXi& T, bool only_interior )
	{
		if( only_interior )
		{
			const auto skip_tet = [ &mesh ]( const int i ) { return mesh.tets[ i ].is_outside; };
			const auto skip_vertex = [ &mesh ]( const int i ) { return mesh.tet_vertices[ i ].is_outside; };
			floatTetWild::extract_volume_mesh( mesh, skip_tet, skip_vertex, V, T );
		}
		else
		{
			const auto skip_tet = [ &mesh ]( const int i ) { return mesh.tets[ i ].is_removed; };
			const auto skip_vertex = [ &mesh ]( const int i ) { return mesh.tet_vertices[ i ].is_removed; };
			floatTetWild::extract_volume_mesh( mesh, skip_tet, skip_vertex, V, T );
		}
	}
}  // namespace floatTetWild
