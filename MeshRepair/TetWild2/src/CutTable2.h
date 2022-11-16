#pragma once
#include <stdint.h>
#include "Types.hpp"

namespace floatTetWild
{
	struct CutTable2
	{
		template<class Store, class View, class Impl>
		class BufferViewBase
		{
			const Store* const m_begin;
			const Store* const m_end;

		  public:
			BufferViewBase( const Store* p1, const Store* p2 )
				: m_begin( p1 )
				, m_end( p2 )
			{
			}

			size_t size() const
			{
				return m_end - m_begin;
			}
			bool empty() const
			{
				return m_end <= m_begin;
			}

			class Iterator
			{
				const Store* pointer;

			  public:
				Iterator( const Store* p )
					: pointer( p )
				{
				}
				void operator++()
				{
					pointer++;
				}
				bool operator!=( Iterator b ) const
				{
					return pointer != b.pointer;
				}
				View operator*() const
				{
					return Impl::loadElement( pointer );
				}
				const Store* getPointer() const
				{
					return pointer;
				}
			};

			Iterator begin() const
			{
				return Iterator( m_begin );
			}
			Iterator end() const
			{
				return Iterator( m_end );
			}

			View operator[]( size_t idx ) const
			{
				assert( idx < size() );
				return Impl::loadElement( m_begin + idx );
			}
		};

		class Vec4Buffer : public BufferViewBase<uint32_t, Vector4i, Vec4Buffer>
		{
		  public:
			using BufferViewBase::BufferViewBase;

			static inline Vector4i loadElement( const uint32_t* rsi )
			{
				// Following 2 lines compile into a single instruction, because pmovsxbd supports memory operand:
				// https://www.felixcloutier.com/x86/pmovsx
				__m128i v = _mm_loadu_si32( rsi );
				v = _mm_cvtepi8_epi32( v );

				// Assuming the compiler gonna inline things, trying to avoid stores into a temporary object
				return Vector4i {
				  _mm_cvtsi128_si32( v ),
				  _mm_extract_epi32( v, 1 ),
				  _mm_extract_epi32( v, 2 ),
				  _mm_extract_epi32( v, 3 ),
				};
			}
		};

		class Vec2Buffer : public BufferViewBase<uint16_t, Vector2i, Vec2Buffer>
		{
		  public:
			using BufferViewBase::BufferViewBase;

			static inline Vector2i loadElement( const uint16_t* rsi )
			{
				const uint8_t* p = (const uint8_t*)rsi;
				return Vector2i { p[ 0 ], p[ 1 ] };
			}

			bool equal( const Vector2i* beginPtr, const Vector2i* endPtr ) const;
		};

		class Bool4Buffer : public BufferViewBase<uint8_t, std::array<bool, 4>, Bool4Buffer>
		{
		  public:
			using BufferViewBase::BufferViewBase;

			static inline std::array<bool, 4> loadElement( const uint8_t* rsi )
			{
				const uint8_t b = *rsi;
				return std::array<bool, 4> {
				  0 != ( b & 1 ),
				  0 != ( b & 2 ),
				  0 != ( b & 4 ),
				  0 != ( b & 8 ),
				};
			}
		};

		class TetConfigs : public BufferViewBase<uint16_t, Vec4Buffer, TetConfigs>
		{
		  public:
			using BufferViewBase::BufferViewBase;
			static Vec4Buffer loadElement( const uint16_t* rsi );
		};
		class TetConfigsOuter;

		static const TetConfigs get_tet_confs( const int idx );

		class DiagConfigs : public BufferViewBase<uint16_t, Vec2Buffer, DiagConfigs>
		{
		  public:
			using BufferViewBase::BufferViewBase;
			static Vec2Buffer loadElement( const uint16_t* rsi );
		};
		class DiagConfigsOuter;

		// The maximum size of these buffers is 4
		static const DiagConfigs get_diag_confs( const int idx );

		class SurfaceConfigs : public BufferViewBase<uint16_t, Bool4Buffer, SurfaceConfigs>
		{
		  public:
			using BufferViewBase::BufferViewBase;
			static Bool4Buffer loadElement( const uint16_t* rsi );
		};
		class SurfaceConfigsOuter;

		static const SurfaceConfigs get_surface_conf( const int idx );

		class FaceIdConfigs : public BufferViewBase<uint16_t, Vec4Buffer, FaceIdConfigs>
		{
		  public:
			using BufferViewBase::BufferViewBase;
			static Vec4Buffer loadElement( const uint16_t* rsi );
		};
		class FaceIdConfigsOuter;
		static const FaceIdConfigs get_face_id_conf( const int idx );

		// The maximum size of these buffers is 16
		static inline const Vec4Buffer get_tet_conf( const int idx, const int cfg )
		{
			return get_tet_confs( idx )[ cfg ];
		}
		static inline const Bool4Buffer get_surface_conf( const int idx, const int cfg )
		{
			return get_surface_conf( idx )[ cfg ];
		}
		static inline const Vec4Buffer get_face_id_conf( const int idx, const int cfg )
		{
			return get_face_id_conf( idx )[ cfg ];
		}
	};

	using CutTable = CutTable2;
}  // namespace floatTetWild