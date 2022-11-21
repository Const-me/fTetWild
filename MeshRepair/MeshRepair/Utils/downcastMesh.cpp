#include "stdafx.h"
#include "downcastMesh.h"
#include "../../TetWild2/Utils/AvxMath.h"
#include <unordered_map>

namespace
{
	struct VertexKey
	{
		std::array<uint32_t, 3> val;
		VertexKey()
		{
		}
		VertexKey( __m128 pos )
		{
			_mm_store_sd( (double*)&val[ 0 ], _mm_castps_pd( pos ) );
			val[ 2 ] = (uint32_t)_mm_extract_ps( pos, 2 );
		}

		bool operator==( const VertexKey& that ) const
		{
			uint64_t xy1 = *(const uint64_t*)val.data();
			uint64_t xy2 = *(const uint64_t*)that.val.data();
			uint64_t z1 = val[ 2 ];
			uint64_t z2 = that.val[ 2 ];
			return 0 == ( ( xy1 ^ xy2 ) | ( z1 ^ z2 ) );
		}
	};

	struct VertexHash
	{
		size_t operator()( const VertexKey& pos ) const
		{
			uint64_t xy = *(const uint64_t*)pos.val.data();
			return ( 17 * 31 * 31 ) + ( xy * 31 ) + pos.val[ 2 ];
		}
	};

	__forceinline __m128 removeNegativeZeros( __m128 v )
	{
		const __m128 zz = _mm_cmpeq_ps( v, _mm_setzero_ps() );
		return _mm_andnot_ps( zz, v );
	}

	class VertexUndupe
	{
		std::unordered_map<VertexKey, uint32_t, VertexHash> hash;

	  public:
		__forceinline uint32_t addVertex( __m256d v )
		{
			__m128 pos = _mm256_cvtpd_ps( v );
			pos = removeNegativeZeros( pos );
			uint32_t newId = (uint32_t)hash.size();

			auto res = hash.try_emplace( pos, newId );
			return res.first->second;
		}

		size_t size() const
		{
			return hash.size();
		}

		void serialize( std::vector<std::array<float, 3>>& result ) const
		{
			result.resize( hash.size() );
			VertexKey* const rdi = (VertexKey*)result.data();
			for( const auto& p : hash )
				rdi[ p.second ] = p.first;
		}

		void clear()
		{
			std::unordered_map<VertexKey, uint32_t, VertexHash> hh;
			hash.swap( hh );
		}
	};
}  // namespace

void downcastFloats( float* rdi, size_t count, const double* rsi )
{
	const double* const rsiEnd = rsi + count;
	const double* const rsiEndAligned = rsi + ( count / 4 ) * 4;

	for( ; rsi < rsiEndAligned; rsi += 4, rdi += 4 )
	{
		__m256d vd = _mm256_loadu_pd( rsi );
		__m128 vf = _mm256_cvtpd_ps( vd );
		_mm_storeu_ps( rdi, vf );
	}

#pragma loop( no_vector )
	for( ; rsi < rsiEnd; rsi++, rdi++ )
	{
		double d = *rsi;
		float f = (float)d;
		*rdi = f;
	}
}

void downcastMesh( const floatTetWild::VertexBuffer& vertices, const floatTetWild::SurfaceIndexBuffer& faces, std::vector<std::array<float, 3>>& vb32,
  std::vector<std::array<uint32_t, 3>>& ib32 )
{
	const size_t countVerts = vertices.rows();
	const double* rsi = vertices.data();
	const double* const rsiEnd = rsi + countVerts * 3;
	const double* const rsiEndAligned = rsiEnd - 3;
	VertexUndupe undupe;
	std::vector<uint32_t> remap;
	remap.resize( countVerts, UINT_MAX );

	size_t vi = 0;
	for( ; rsi < rsiEndAligned; rsi += 3, vi++ )
	{
		__m256d v = _mm256_loadu_pd( rsi );
		remap[ vi ] = undupe.addVertex( v );
	}

	// For the last vertex we can't use full-vector unaligned load, may fail in runtime with access violation
	{
		__m256d v = AvxMath::loadDouble3( rsi );
		remap[ vi ] = undupe.addVertex( v );
	}

	const size_t countTriangles = faces.rows();
	if( undupe.size() == countVerts )
	{
		// Downcasting didn't produce duplicate vertices in the mesh
		undupe.clear();
		remap.clear();
		remap.shrink_to_fit();

		// Downcast the whole vertex buffer with AVX instructions
		vb32.resize( countVerts );
		downcastFloats( (float*)vb32.data(), countVerts * 3, vertices.data() );

		// The mesh topology has not changed, memcpy will do
		ib32.resize( countTriangles );
		memcpy( ib32.data(), faces.data(), countTriangles * 12 );

		return;
	}

	undupe.serialize( vb32 );
	undupe.clear();

	ib32.resize( countTriangles );
	std::array<uint32_t, 3>* rdi = ib32.data();
	const uint32_t* pFace = (const uint32_t*)faces.data();
	const uint32_t* const pFaceEnd = pFace + countTriangles * 3;
	while( pFace < pFaceEnd )
	{
		uint32_t i0 = pFace[ 0 ];
		uint32_t i1 = pFace[ 1 ];
		uint32_t i2 = pFace[ 2 ];
		pFace += 3;

		i0 = remap[ i0 ];
		i1 = remap[ i1 ];
		i2 = remap[ i2 ];

		if( i0 != i1 && i0 != i2 && i1 != i2 )
		{
			std::array<uint32_t, 3>& res = *rdi;
			res[ 0 ] = i0;
			res[ 1 ] = i1;
			res[ 2 ] = i2;
			rdi++;
		}
	}

	ib32.resize( rdi - ib32.data() );
}