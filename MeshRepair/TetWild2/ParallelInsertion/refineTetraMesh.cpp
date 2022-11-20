#include "stdafx.h"
#include "refineTetraMesh.h"
#include "../Utils/AvxMath.h"
#include <unordered_set>

namespace
{
	// We don't need precision and we want to minimize the count of extra vertices in the mesh
	static const double tetraBoxEpsilonMul = 257.0 / 256.0;
	constexpr size_t maxPasses = 5;

	using namespace floatTetWild;

	struct RefineContext
	{
		std::vector<Vector3>& vertices;
		iDelaunay& delaunay;

		void initialize( const Vector3& boxMin, const Vector3& boxMax, __m128i voxels, const std::vector<Vector3i>& inputFaces );

		__m256d refine( bool insertVertices, bool* insertedSome = nullptr ) const;

		RefineContext( std::vector<Vector3>& vb, iDelaunay& del )
			: vertices( vb )
			, delaunay( del )
		{
		}

	  private:
		double lengthThresholdSq = 0;
		std::vector<uint64_t> sourceEdges;
		mutable std::unordered_set<uint64_t> splitEdges;

		inline bool isSourceEdge( uint64_t key ) const;

		template<int i0, int i1>
		inline void trySplitEdge( double len, double th, __m128i element, __m256d v0, __m256d v1, bool& original ) const;

		void trySplitEdge( double len, double th, uint32_t i0, uint32_t i1, __m256d v0, __m256d v1, bool& original ) const;
	};

	__forceinline uint64_t edgeKey( uint64_t a, uint64_t b )
	{
		uint64_t u = ( a << 32 ) | b;
		uint64_t flip = _rotr64( u, 32 );
		return std::max( u, flip );
	}

	static void collectEdgeIDs( std::vector<uint64_t>& res, const std::vector<Vector3i>& faces )
	{
		res.resize( faces.size() * 3 );
		uint64_t* rdi = res.data();
		for( const Vector3i& face : faces )
		{
			const uint64_t a = (uint32_t)face[ 0 ];
			const uint64_t b = (uint32_t)face[ 1 ];
			const uint64_t c = (uint32_t)face[ 2 ];
			rdi[ 0 ] = edgeKey( a, b );
			rdi[ 1 ] = edgeKey( a, c );
			rdi[ 2 ] = edgeKey( b, c );
			rdi += 3;
		}

		std::sort( res.begin(), res.end() );
		res.erase( std::unique( res.begin(), res.end() ), res.end() );
		res.shrink_to_fit();
	}

	void RefineContext::initialize( const Vector3& boxMin, const Vector3& boxMax, __m128i voxels, const std::vector<Vector3i>& inputFaces )
	{
		using namespace AvxMath;
		__m256d min = loadDouble3( boxMin.data() );
		__m256d max = loadDouble3( boxMax.data() );
		__m256d size = _mm256_sub_pd( max, min );

		voxels = _mm_shuffle_epi32( voxels, _MM_SHUFFLE( 2, 2, 1, 0 ) );
		__m256d div = _mm256_cvtepi32_pd( voxels );
		__m256d res = _mm256_div_pd( size, div );

		__m256d maxBc = vector3BroadcastMaximum( res );
		double t = _mm256_cvtsd_f64( maxBc );
		t *= tetraBoxEpsilonMul;
		lengthThresholdSq = t * t;

		collectEdgeIDs( sourceEdges, inputFaces );
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

	inline double edgeLength( __m256d v0, __m256d v1, double threshold, uint32_t& bmp, uint32_t bit )
	{
		__m256d dist = _mm256_sub_pd( v1, v0 );
		double lsq = AvxMath::vector3DotScalar( dist, dist );
		bmp |= ( lsq > threshold ) ? bit : 0;
		return lsq;
	}

	template<int i>
	inline uint32_t extractLane( __m128i v )
	{
		static_assert( i >= 0 && i < 4 );
		if constexpr( i == 0 )
			return (uint32_t)_mm_cvtsi128_si32( v );
		else
			return (uint32_t)_mm_extract_epi32( v, i );
	}

	inline bool RefineContext::isSourceEdge( uint64_t key ) const
	{
		// That vector is sorted, can use binary search instead of linear
		auto it = std::lower_bound( sourceEdges.begin(), sourceEdges.end(), key );
		return it != sourceEdges.end() && *it == key;
	}

	void RefineContext::trySplitEdge( double len, double th, uint32_t i0, uint32_t i1, __m256d v0, __m256d v1, bool& original ) const
	{
		const uint64_t key = edgeKey( i0, i1 );

		if( isSourceEdge( key ) )
		{
			// Avoid splitting edges of the original mesh, introduces new vertices very close to the original surface
			// For that case, we instead inserting new vertex inside the element, in the center of mass
			original = true;
			return;
		}

		if( !splitEdges.emplace( key ).second )
		{
			// That edge is shared by another element, and we already split that very edge while processing another element on the current iteration
			return;
		}

		__m128d tmp = _mm_set_sd( len / th );	 // length^2 / threshold^2
		tmp = _mm_sqrt_sd( tmp, tmp );			 // take a square root
		tmp = _mm_ceil_sd( tmp, tmp );			 // Round up to an integer, SSE 4.1 has a special instruction for that
		const int count = _mm_cvtsd_i32( tmp );	 // Probably round to nearest, because at least on Windows that's in the default state of MXCSR register
		if( count < 2 )
			return;	 // Probably FP precision shenanigans?

		// Produce midpoint[s] on the edge being split
		const double div = (double)count;
		for( int i = 1; i < count; i++ )
		{
			const __m256d pos = AvxMath::lerpFast( v0, v1, (double)i / div );
			Vector3& dest = vertices.emplace_back();
			AvxMath::storeDouble3( dest.data(), pos );
		}
	}

	template<int i0, int i1>
	inline void RefineContext::trySplitEdge( double len, double th, __m128i element, __m256d v0, __m256d v1, bool& original ) const
	{
		static_assert( i0 != i1 );
		uint32_t u0 = extractLane<i0>( element );
		uint32_t u1 = extractLane<i1>( element );
		trySplitEdge( len, th, u0, u1, v0, v1, original );
	}

	__m256d RefineContext::refine( bool insertVertices, bool* insertedSome ) const
	{
		splitEdges.clear();
		bool anyNewVertex = false;
		const __m128i* rsi = delaunay.getElements();
		const __m128i* const rsiEnd = rsi + delaunay.countElements();
		__m256d res = _mm256_setzero_pd();

		alignas( 32 ) std::array<double, 6> edgeLengths;
		const double th = lengthThresholdSq;
		for( ; rsi < rsiEnd; rsi++ )
		{
			const __m128i element = *rsi;
			__m256d v0, v1, v2, v3;
			loadElement( vertices, element, v0, v1, v2, v3 );
			const __m256d es = sizeElement( v0, v1, v2, v3 );
			res = _mm256_max_pd( res, es );
			if( !insertVertices )
				continue;

			uint32_t bitmap = 0;
			edgeLengths[ 0 ] = edgeLength( v0, v1, th, bitmap, 1 );
			edgeLengths[ 1 ] = edgeLength( v0, v2, th, bitmap, 2 );
			edgeLengths[ 2 ] = edgeLength( v0, v3, th, bitmap, 4 );
			edgeLengths[ 3 ] = edgeLength( v1, v2, th, bitmap, 8 );
			edgeLengths[ 4 ] = edgeLength( v1, v3, th, bitmap, 0x10 );
			edgeLengths[ 5 ] = edgeLength( v2, v3, th, bitmap, 0x20 );
			if( 0 == bitmap )
				continue;  // All edges are within the threshold

			anyNewVertex = true;

			bool anySourceEdge = false;
			if( 0 != ( bitmap & 1 ) )
				trySplitEdge<0, 1>( edgeLengths[ 0 ], th, element, v0, v1, anySourceEdge );
			if( 0 != ( bitmap & 2 ) )
				trySplitEdge<0, 2>( edgeLengths[ 1 ], th, element, v0, v2, anySourceEdge );
			if( 0 != ( bitmap & 4 ) )
				trySplitEdge<0, 3>( edgeLengths[ 2 ], th, element, v0, v3, anySourceEdge );
			if( 0 != ( bitmap & 8 ) )
				trySplitEdge<1, 2>( edgeLengths[ 3 ], th, element, v1, v2, anySourceEdge );
			if( 0 != ( bitmap & 0x10 ) )
				trySplitEdge<1, 3>( edgeLengths[ 4 ], th, element, v1, v3, anySourceEdge );
			if( 0 != ( bitmap & 0x20 ) )
				trySplitEdge<2, 3>( edgeLengths[ 5 ], th, element, v2, v3, anySourceEdge );

			if( anySourceEdge )
			{
				const __m256d center = elementCenter( v0, v1, v2, v3 );
				Vector3& nv = vertices.emplace_back();
				AvxMath::storeDouble3( nv.data(), center );
			}
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

__m256d floatTetWild::refineTetraMesh(
  const Vector3& boxMin, const Vector3& boxMax, iDelaunay& delaunay, std::vector<Vector3>& vertices, __m128i voxels, const std::vector<Vector3i>& inputFaces )
{
	RefineContext context { vertices, delaunay };

	context.initialize( boxMin, boxMax, voxels, inputFaces );

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