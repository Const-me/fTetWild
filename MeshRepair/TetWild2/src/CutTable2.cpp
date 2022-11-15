#include "stdafx.h"
#include "CutTable2.h"

namespace
{
	namespace Data
	{
#include "CutTable2.inl"
	}
}  // namespace

namespace floatTetWild
{
	CutTable2::Vec4Buffer CutTable2::TetConfigs::loadElement( const uint16_t* rsi )
	{
		const uint32_t* const payload = Data::TetraConfigsData.data();
		return Vec4Buffer {
		  payload + rsi[ 0 ],
		  payload + rsi[ 1 ],
		};
	}
	class CutTable2::TetConfigsOuter : public BufferViewBase<uint8_t, TetConfigs, TetConfigsOuter>
	{
	  public:
		using BufferViewBase::BufferViewBase;
		static TetConfigs loadElement( const uint8_t* rsi )
		{
			const uint16_t* const payload = Data::TetraConfigsInner.data();
			return TetConfigs {
			  payload + rsi[ 0 ],
			  payload + rsi[ 1 ],
			};
		}
	};
	static const CutTable2::TetConfigsOuter g_tetConfigs {
	  Data::TetraConfigsOuter.data(),
	  Data::TetraConfigsOuter.data() + ( Data::TetraConfigsOuter.size() - 1 ),
	};
	const CutTable2::TetConfigs CutTable2::get_tet_confs( const int idx )
	{
		return g_tetConfigs[ (size_t)idx ];
	}

	CutTable2::Vec2Buffer CutTable2::DiagConfigs::loadElement( const uint16_t* rsi )
	{
		const uint16_t* const payload = Data::DiagConfigsData.data();
		return Vec2Buffer {
		  payload + rsi[ 0 ],
		  payload + rsi[ 1 ],
		};
	}
	class CutTable2::DiagConfigsOuter : public BufferViewBase<uint8_t, DiagConfigs, DiagConfigsOuter>
	{
	  public:
		using BufferViewBase::BufferViewBase;
		static DiagConfigs loadElement( const uint8_t* rsi )
		{
			const uint16_t* const payload = Data::DiagConfigsInner.data();
			return DiagConfigs {
			  payload + rsi[ 0 ],
			  payload + rsi[ 1 ],
			};
		}
	};
	static const CutTable2::DiagConfigsOuter g_diagConfigs {
	  Data::DiagConfigsOuter.data(),
	  Data::DiagConfigsOuter.data() + Data::DiagConfigsOuter.size() - 1,
	};
	const CutTable2::DiagConfigs CutTable2::get_diag_confs( const int idx )
	{
		return g_diagConfigs[ (size_t)idx ];
	}

	CutTable2::Bool4Buffer CutTable2::SurfaceConfigs::loadElement( const uint16_t* rsi )
	{
		const uint8_t* const payload = Data::SurfaceConfigData.data();
		return Bool4Buffer {
		  payload + rsi[ 0 ],
		  payload + rsi[ 1 ],
		};
	}
	class CutTable2::SurfaceConfigsOuter : public BufferViewBase<uint8_t, SurfaceConfigs, SurfaceConfigsOuter>
	{
	  public:
		using BufferViewBase::BufferViewBase;
		static SurfaceConfigs loadElement( const uint8_t* rsi )
		{
			const uint16_t* const payload = Data::SurfaceConfigInner.data();
			return SurfaceConfigs {
			  payload + rsi[ 0 ],
			  payload + rsi[ 1 ],
			};
		}
	};
	static const CutTable2::SurfaceConfigsOuter g_surfaceConfigs {
	  Data::SurfaceConfigOuter.data(),
	  Data::SurfaceConfigOuter.data() + Data::SurfaceConfigOuter.size() - 1,
	};
	const CutTable2::SurfaceConfigs CutTable2::get_surface_conf( const int idx )
	{
		return g_surfaceConfigs[ (size_t)idx ];
	}

	CutTable2::Vec4Buffer CutTable2::FaceIdConfigs::loadElement( const uint16_t* rsi )
	{
		const uint32_t* const payload = Data::FaceIDsData.data();
		return Vec4Buffer {
		  payload + rsi[ 0 ],
		  payload + rsi[ 1 ],
		};
	}
	class CutTable2::FaceIdConfigsOuter : public BufferViewBase<uint8_t, FaceIdConfigs, FaceIdConfigsOuter>
	{
	  public:
		using BufferViewBase::BufferViewBase;
		static FaceIdConfigs loadElement( const uint8_t* rsi )
		{
			const uint16_t* const payload = Data::FaceIDsInner.data();
			return FaceIdConfigs {
			  payload + rsi[ 0 ],
			  payload + rsi[ 1 ],
			};
		}
	};
	static const CutTable2::FaceIdConfigsOuter g_faceIdConfigs {
	  Data::FaceIDsOuter.data(),
	  Data::FaceIDsOuter.data() + Data::FaceIDsOuter.size() - 1,
	};
	const CutTable2::FaceIdConfigs CutTable2::get_face_id_conf( const int idx )
	{
		return g_faceIdConfigs[ (size_t)idx ];
	}

	bool CutTable2::Vec2Buffer::equal( const Vector2i* beginPtr, const Vector2i* endPtr ) const
	{
		const size_t length = endPtr - beginPtr;
		if( length != size() )
			return false;
		if( 0 == length )
			return true;

		const uint16_t* p1 = begin().getPointer();
		const Vector2i* p2 = beginPtr;
		const uint16_t* const p1End = p1 + length;
		const uint16_t* const p1EndAligned = p1 + ( length / 2 ) * 2;

		for( ; p1 < p1EndAligned; p1 += 2, p2 += 2 )
		{
			__m128i v1 = _mm_loadu_si32( p1 );
			v1 = _mm_cvtepu8_epi32( v1 );
			__m128i v2 = _mm_loadu_si128( (const __m128i*)p2 );
			__m128i xx = _mm_xor_si128( v1, v2 );
			if( _mm_testz_si128( xx, xx ) )
				continue;
			return false;
		}

		if( p1 < p1End )
		{
			Vector2i e1 = CutTable2::Vec2Buffer::loadElement( p1 );
			if( e1 != *p2 )
				return false;
		}
		return true;
	}
}  // namespace floatTetWild