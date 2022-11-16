// This file was originally part of libigl, a simple c++ geometry processing library,
// subject to the terms of the Mozilla Public License v. 2.0.
// I've replaced the implementation with another one, based on the OpenMP thread pool.
#pragma once
#include "igl_inline.h"
#include <functional>
#include "parallelThreadsCount.h"

namespace igl
{

	// PARALLEL_FOR Functional implementation of a basic, open-mp style, parallel
	// for loop. If the inner block of a for-loop can be rewritten/encapsulated in
	// a single (anonymous/lambda) function call `func` so that the serial code
	// looks like:
	//
	//     for(int i = 0;i<loop_size;i++)
	//     {
	//       func(i);
	//     }
	//
	// then `parallel_for(loop_size,func,min_parallel)` will use as many threads as
	// available on the current hardware to parallelize this for loop so long as
	// loop_size<min_parallel, otherwise it will just use a serial for loop.
	//
	// Inputs:
	//   loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
	//   func  function handle taking iteration index as only argument to compute
	//     inner block of for loop I.e. for(int i ...){ func(i); }
	//   min_parallel  min size of loop_size such that parallel (non-serial)
	//     thread pooling should be attempted {0}
	// Returns true iff thread pool was invoked
	template<typename Index, typename FunctionType>
	inline bool parallel_for( const Index loop_size, const FunctionType& func, const size_t min_parallel = 0 );
	// PARALLEL_FOR Functional implementation of an open-mp style, parallel for
	// loop with accumulation. For example, serial code separated into n chunks
	// (each to be parallelized with a thread) might look like:
	//
	//     Eigen::VectorXd S;
	//     const auto & prep_func = [&S](int n){ S = Eigen:VectorXd::Zero(n); };
	//     const auto & func = [&X,&S](int i, int t){ S(t) += X(i); };
	//     const auto & accum_func = [&S,&sum](int t){ sum += S(t); };
	//     prep_func(n);
	//     for(int i = 0;i<loop_size;i++)
	//     {
	//       func(i,i%n);
	//     }
	//     double sum = 0;
	//     for(int t = 0;t<n;t++)
	//     {
	//       accum_func(t);
	//     }
	//
	// Inputs:
	//   loop_size  number of iterations. I.e. for(int i = 0;i<loop_size;i++) ...
	//   prep_func function handle taking n >= number of threads as only
	//     argument
	//   func  function handle taking iteration index i and thread id t as only
	//     arguments to compute inner block of for loop I.e.
	//     for(int i ...){ func(i,t); }
	//   accum_func  function handle taking thread index as only argument, to be
	//     called after all calls of func, e.g., for serial accumulation across
	//     all n (potential) threads, see n in description of prep_func.
	//   min_parallel  min size of loop_size such that parallel (non-serial)
	//     thread pooling should be attempted {0}
	// Returns true iff thread pool was invoked
	template<typename Index, typename PrepFunctionType, typename FunctionType, typename AccumFunctionType>
	inline bool parallel_for(
	  const Index loop_size, const PrepFunctionType& prep_func, const FunctionType& func, const AccumFunctionType& accum_func, const size_t min_parallel = 0 );
}  // namespace igl

// Implementation
#include <omp.h>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

template<typename Index, typename FunctionType>
inline bool igl::parallel_for( const Index loop_size, const FunctionType& func, const size_t min_parallel )
{
	using namespace std;
	// no op preparation/accumulation
	const auto& no_op = []( const size_t /*n/t*/ ) {};
	// two-parameter wrapper ignoring thread id
	const auto& wrapper = [ &func ]( Index i, size_t /*t*/ ) { func( i ); };
	return parallel_for( loop_size, no_op, wrapper, no_op, min_parallel );
}

template<typename Index, typename PreFunctionType, typename FunctionType, typename AccumFunctionType>
inline bool igl::parallel_for(
  const Index loop_size, const PreFunctionType& prep_func, const FunctionType& func, const AccumFunctionType& accum_func, const size_t min_parallel )
{
	assert( loop_size >= 0 );
	if( loop_size == 0 )
		return false;
		// Estimate number of threads in the pool
		// http://ideone.com/Z7zldb
#ifdef IGL_PARALLEL_FOR_FORCE_SERIAL
	const size_t nthreads = 1;
#else
	const size_t nthreads = MeshRepair::getThreadsCount();
#endif
	if( loop_size < min_parallel || nthreads <= 1 )
	{
		// serial
		prep_func( 1 );
		for( Index i = 0; i < loop_size; i++ )
			func( i, 0 );
		accum_func( 0 );
		return false;
	}
	else
	{
		// Size of a slice for the range functions
		Index slice = std::max( (Index)std::round( ( loop_size + 1 ) / static_cast<double>( nthreads ) ), (Index)1 );

		// [Helper] Inner loop
		const auto& range = [ &func ]( const Index k1, const Index k2, const size_t t )
		{
			for( Index k = k1; k < k2; k++ )
				func( k, t );
		};
		prep_func( nthreads );
		// Create pool and launch jobs
		struct ThreadJob
		{
			Index k1, k2;
			size_t t;
		};
		std::vector<ThreadJob> pool;
		pool.reserve( nthreads );
		// Inner range extents
		Index i1 = 0;
		Index i2 = std::min( 0 + slice, loop_size );
		{
			size_t t = 0;
			for( ; t + 1 < nthreads && i1 < loop_size; ++t )
			{
				pool.emplace_back( ThreadJob { i1, i2, t } );
				i1 = i2;
				i2 = std::min( i2 + slice, loop_size );
			}
			if( i1 < loop_size )
				pool.emplace_back( ThreadJob { i1, loop_size, t } );
		}

#pragma omp parallel for schedule( dynamic, 1 ) 
		for( int i = 0; i < (int)pool.size(); i++ )
		{
			const ThreadJob& j = pool[ i ];
			range( j.k1, j.k2, j.t );
		}

		// Accumulate across threads
		for( size_t t = 0; t < nthreads; t++ )
			accum_func( t );
		return true;
	}
}