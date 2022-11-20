#include "stdafx.h"
#include "refineTetraMesh.h"
#include "../Utils/AvxMath.h"

namespace
{
	// We don't need precision and we want to minimize the count of extra vertices in the mesh
	static const double tetraBoxEpsilonMul = 257.0 / 256.0;
	constexpr size_t maxPasses = 5;

	using namespace floatTetWild;

	struct RefineContext
	{
		const __m256d threshold;
		std::vector<Vector3>& vertices;
		iDelaunay& delaunay;

		__m256d refine( bool insertVertices, bool* insertedSome = nullptr ) const;
	};

	__m256d tetraSizeThreshold( const Vector3& boxMin, const Vector3& boxMax, __m128i voxels )
	{
		using namespace AvxMath;
		__m256d min = loadDouble3( boxMin.data() );
		__m256d max = loadDouble3( boxMax.data() );
		__m256d size = _mm256_sub_pd( max, min );

		voxels = _mm_shuffle_epi32( voxels, _MM_SHUFFLE( 2, 2, 1, 0 ) );
		__m256d div = _mm256_cvtepi32_pd( voxels );
		__m256d res = _mm256_div_pd( size, div );

		__m256d safetyMul = _mm256_broadcast_sd( &tetraBoxEpsilonMul );
		return _mm256_mul_pd( res, safetyMul );
	}

	inline __m256d loadVertex( const std::vector<Vector3>& vertices, int si )
	{
		const size_t ui = (size_t)(uint32_t)si;
		const double* const p = vertices[ ui ].data();
		if( ui + 1 < vertices.size() )
			return _mm256_loadu_pd( p );  //< This branch is extremely likely to be taken, and saves quite a few instructions
		else
			return AvxMath::loadDouble3( p );
	}

	inline void loadElement( const std::vector<Vector3>& vertices, __m128i elt, __m256d& v0, __m256d& v1, __m256d& v2, __m256d& v3 )
	{
		v0 = loadVertex( vertices, _mm_cvtsi128_si32( elt ) );
		v1 = loadVertex( vertices, _mm_extract_epi32( elt, 1 ) );
		v2 = loadVertex( vertices, _mm_extract_epi32( elt, 2 ) );
		v3 = loadVertex( vertices, _mm_extract_epi32( elt, 3 ) );
	}

	inline __m256d sizeElement( __m256d v0, __m256d v1, __m256d v2, __m256d v3 )
	{
		__m256d i = _mm256_min_pd( _mm256_min_pd( v0, v1 ), _mm256_min_pd( v2, v3 ) );
		__m256d ax = _mm256_max_pd( _mm256_max_pd( v0, v1 ), _mm256_max_pd( v2, v3 ) );
		return _mm256_sub_pd( ax, i );
	}

	inline __m256d elementCenter( __m256d v0, __m256d v1, __m256d v2, __m256d v3 )
	{
		__m256d a = _mm256_add_pd( v0, v1 );
		__m256d b = _mm256_add_pd( v2, v3 );
		__m256d r = _mm256_add_pd( a, b );
		return _mm256_mul_pd( r, _mm256_set1_pd( 0.25 ) );
	}

	__m256d RefineContext::refine( bool insertVertices, bool* insertedSome ) const
	{
		bool anyNewVertex = false;
		const __m128i* rsi = delaunay.getElements();
		const __m128i* const rsiEnd = rsi + delaunay.countElements();
		__m256d res = _mm256_setzero_pd();
		for( ; rsi < rsiEnd; rsi++ )
		{
			__m256d v0, v1, v2, v3;
			loadElement( vertices, *rsi, v0, v1, v2, v3 );
			const __m256d es = sizeElement( v0, v1, v2, v3 );
			res = _mm256_max_pd( res, es );
			if( !insertVertices )
				continue;

			const __m256d gt = _mm256_cmp_pd( es, threshold, _CMP_GT_OQ );
			uint32_t mask = (uint32_t)_mm256_movemask_pd( gt );
			if( 0 == ( mask & 0b111 ) )
				continue;

			const __m256d center = elementCenter( v0, v1, v2, v3 );
			Vector3& nv = vertices.emplace_back();
			AvxMath::storeDouble3( nv.data(), center );
			anyNewVertex = true;
		}

		if( nullptr != insertedSome )
			*insertedSome = anyNewVertex;
		return res;
	}

	static const double one = 1;
	inline __m256d setOneW( __m256d xyz )
	{
		__m256d w = _mm256_broadcast_sd( &one );
		return _mm256_blend_pd( xyz, w, 0b1000 );
	}
}  // namespace

__m256d floatTetWild::refineTetraMesh( const Vector3& boxMin, const Vector3& boxMax, iDelaunay& delaunay, std::vector<Vector3>& vertices, __m128i voxels )
{
	RefineContext context { tetraSizeThreshold( boxMin, boxMax, voxels ), vertices, delaunay };

	for( size_t pass = 0; pass < maxPasses; pass++ )
	{
		bool inserted;
		__m256d res = context.refine( true, &inserted );
		if( !inserted )
			return setOneW( res );
		delaunay.compute( vertices.size(), (const double*)vertices.data() );
	}

	return setOneW( context.refine( false, nullptr ) );
}