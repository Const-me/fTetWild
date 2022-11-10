#include "stdafx.h"
#include "meshRepairMain.h"
#include "convertParameters.h"
#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif
#include "../TetWild2/src/AABBWrapper.h"
#include "../TetWild2/src/Simplification.h"
#include "../TetWild2/src/FloatTetDelaunay.h"
#include "../TetWild2/src/TriangleInsertion.h"
#include "../TetWild2/src/MeshImprovement.h"
#include "ResultMesh.h"
#include "Utils/writeStl.h"

using namespace floatTetWild;

namespace
{
	constexpr bool skip_simplify = false;
}

HRESULT meshRepairMain(
  MeshRepair::SourceMesh& rsi, const MeshRepair::Parameters& parameters, const MeshRepair::sLoggerSetup& logger, MeshRepair::iResultMesh** rdi ) noexcept
{
	// writeStl( rsi.input_vertices, rsi.input_faces, LR"(C:\Temp\2remove\MeshRepair\Temp-01.stl)" );

	Mesh mesh { logger };
	try
	{
		Parameters& params = mesh.params;
		CHECK( convertParameters( params, parameters ) );

#ifdef FLOAT_TETWILD_USE_TBB
		uint32_t max_threads = parameters.maxThreads;
		if( 0 == max_threads )
			max_threads = std::thread::hardware_concurrency();

		const size_t MB = 1024 * 1024;
		const size_t stack_size = 64 * MB;
		uint32_t num_threads = std::max( 1u, std::thread::hardware_concurrency() );
		num_threads = std::min( max_threads, num_threads );
		params.num_threads = num_threads;
		std::cout << "TBB threads " << num_threads << std::endl;
		tbb::task_scheduler_init scheduler( num_threads, stack_size );
#endif

		AABBWrapper tree( rsi.mesh );
		if( !params.init( tree.get_sf_diag() ) )
			return E_INVALIDARG;

		// preprocessing
		simplify( rsi.input_vertices, rsi.input_faces, rsi.input_tags, tree, params, skip_simplify );
		tree.init_b_mesh_and_tree( rsi.input_vertices, rsi.input_faces, mesh );

		// tetrahedralizing
		std::vector<bool> is_face_inserted( rsi.input_faces.size(), false );
		FloatTetDelaunay::tetrahedralize( rsi.input_vertices, rsi.input_faces, tree, mesh, is_face_inserted );

		// cutting
		insert_triangles( rsi.input_vertices, rsi.input_faces, rsi.input_tags, mesh, is_face_inserted, tree, false );

		// mesh optimization
		optimization( rsi.input_vertices, rsi.input_faces, rsi.input_tags, is_face_inserted, mesh, tree, { { 1, 1, 1, 1 } } );

		// correct_tracked_surface_orientation
		correct_tracked_surface_orientation( mesh, tree );

		if( params.smooth_open_boundary )
		{
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
				if( params.use_floodfill )
					filter_outside_floodfill( mesh );
				else if( params.use_input_for_wn )
					filter_outside( mesh, rsi.input_vertices, rsi.input_faces );
				else
					filter_outside( mesh );
			}
		}

		using namespace ComLight;
		CComPtr<Object<MeshRepair::ResultMesh>> result;
		CHECK( Object<MeshRepair::ResultMesh>::create( result ) );

		if( params.manifold_surface )
			manifold_surface( mesh, result->V_sf, result->F_sf );
		else
			get_surface( mesh, result->V_sf, result->F_sf );

		result.detach( rdi );
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

#ifdef _MSC_VER
#ifdef FLOAT_TETWILD_USE_TBB
#pragma comment( lib, "tbb_static.lib" )
#endif
#pragma comment( lib, "mpir.lib" )
#endif