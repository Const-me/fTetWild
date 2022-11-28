#include "stdafx.h"
#include "meshRepairMain.h"
#include "convertParameters.h"
#include "../TetWild2/src/AABBWrapper.h"
#include "../TetWild2/src/Simplification.h"
#include "../TetWild2/src/FloatTetDelaunay.h"
#include "../TetWild2/src/TriangleInsertion.h"
#include "../TetWild2/src/MeshImprovement.h"
#include "ResultMesh.h"
// #include "Utils/writeStl.h"
#include "../TetWild2/parallelThreadsImpl.h"
#include "Utils/downcastMesh.h"

using namespace floatTetWild;

HRESULT meshRepairMain( MeshRepair::SourceMesh& rsi, MeshRepair::eGlobalFlags globalFlags, const MeshRepair::Parameters& parameters,
  const MeshRepair::sLoggerSetup& logger, MeshRepair::iResultMesh** rdi ) noexcept
{
	// writeStl( rsi.input_vertices, rsi.input_faces, LR"(C:\Temp\2remove\MeshRepair\Temp-01.stl)" );

	Mesh mesh { logger };
	mesh.logger().logInfo( "Starting mesh repair" );
	using MeshRepair::eRepairFlags;
	try
	{
		Parameters& params = mesh.params;
		CHECK( convertParameters( params, globalFlags, parameters ) );
		MeshRepair::SetThreadsCountRaii iglThreadCountSet { params.num_threads };
		mesh.createThreadLocalBuffers();

		AABBWrapper tree( rsi.mesh, mesh.facetRecursionStacks );
		const double boxDiagonal = tree.get_sf_diag();
		if( parameters.hasFlag( eRepairFlags::LengthsAreAbsolute ) )
		{
			params.ideal_edge_length /= boxDiagonal;
			params.eps_rel /= boxDiagonal;
		}

		params.init( boxDiagonal );

		mesh.logger().logInfo( "Preprocessing.." );
		const bool skipSimplify = parameters.hasFlag( eRepairFlags::SkipSimplify );
		simplify( rsi.input_vertices, rsi.input_faces, rsi.input_tags, tree, params, skipSimplify );
		tree.init_b_mesh_and_tree( rsi.input_vertices, rsi.input_faces, mesh );

		mesh.logger().logInfo( "Creating initial volume mesh.." );
		BoolVector is_face_inserted;
		is_face_inserted.resize( rsi.input_faces.size(), false );
		FloatTetDelaunay::tetrahedralize( rsi.input_vertices, rsi.input_faces, tree, mesh, is_face_inserted );

		mesh.logger().logInfo( "Inserting triangles.." );
		insert_triangles( rsi.input_vertices, rsi.input_faces, rsi.input_tags, mesh, is_face_inserted, tree, false );

		mesh.logger().logInfo( "Optimization passes.." );
		optimization( rsi.input_vertices, rsi.input_faces, rsi.input_tags, is_face_inserted, mesh, tree, { { 1, 1, 1, 1 } } );

		mesh.logger().logInfo( "Correcting tracked surface orientation.." );
		correct_tracked_surface_orientation( mesh, tree );

		if( params.smooth_open_boundary )
		{
			mesh.logger().logInfo( "Smoothing the open boundary.." );
			smooth_open_boundary( mesh, tree );
			for( auto& t : mesh.tets )
			{
				if( t.is_outside )
					t.is_removed = true;
			}
		}
		else
		{
			if( !params.disable_filtering )
			{
				mesh.logger().logInfo( "Filtering out outside elements.." );
				if( params.use_floodfill )
					filter_outside_floodfill( mesh );
				else if( params.use_input_for_wn )
					filter_outside( mesh, rsi.input_vertices, rsi.input_faces );
				else
					filter_outside( mesh );
			}
		}

		VertexBuffer V_sf;
		SurfaceIndexBuffer F_sf;
		if( params.manifold_surface )
		{
			mesh.logger().logInfo( "Generating manifold surface.." );
			manifold_surface( mesh, V_sf, F_sf );
		}
		else
		{
			mesh.logger().logInfo( "Gathering surface mesh.." );
			get_surface( mesh, V_sf, F_sf );
		}

		mesh.times.logDebug( mesh.logger() );

		using namespace ComLight;
		if( 0 == ( parameters.flags & MeshRepair::eRepairFlags::DowncastMeshFp32 ) ) 
		{
			CComPtr<Object<MeshRepair::ResultMesh>> result;
			CHECK( Object<MeshRepair::ResultMesh>::create( result ) );
			result->V_sf = std::move( V_sf );
			result->F_sf = std::move( F_sf );
			result.detach( rdi );
		}
		else 
		{
			mesh.logger().logInfo( "Generating FP32 indexed mesh.." );
			std::vector<std::array<float, 3>> vb32;
			std::vector<std::array<uint32_t, 3>> ib32;
			downcastMesh( V_sf, F_sf, vb32, ib32 );

			CComPtr<Object<MeshRepair::ResultMesh32>> result;
			CHECK( Object<MeshRepair::ResultMesh32>::create( result ) );
			result->vertexBuffer = std::move( vb32 );
			result->indexBuffer = std::move( ib32 );
			result.detach( rdi );
		}
		mesh.logger().logInfo( "Mesh repair complete." );
		return S_OK;
	}
	catch( const std::bad_alloc& )
	{
		return E_OUTOFMEMORY;
	}
	catch( const std::exception& ex )
	{
		mesh.logger().logError( "Mesh repair failed: %s", ex.what() );
		return E_FAIL;
	}
}

// If you try to compile with gcc or clang, you gonna need to add libgmp.so or an equivalent to the input of the linker
#ifdef _MSC_VER
#pragma comment( lib, "mpir.lib" )
#endif