#include "stdafx.h"
#include "cutTableData.h"
#include "../src/auto_table.hpp"

namespace
{
	using floatTetWild::Vector2i;
	using floatTetWild::Vector4i;

	// Do not try this at home, this is a design-time code, quality is irrelevant as long as the result is good
	FILE* g_file;
#define printf( s, ... ) fprintf( g_file, s, __VA_ARGS__ )

	template<class Store, class View>
	class InnerIndexBase
	{
		std::vector<Store> vec;

	  protected:
		virtual const char* storeTypeName() const = 0;
		virtual Store compressElement( const View& elt ) const = 0;
		virtual void printElement( Store elt ) const = 0;

	  public:
		uint32_t size() const
		{
			return (uint32_t)vec.size();
		}

		void add( const View& val )
		{
			vec.push_back( compressElement( val ) );
		}

		void print( const char* name ) const
		{
			constexpr size_t elementsPerLine = 16;

			printf( "static const std::array<%s, %i> %s%s = \n", storeTypeName(), (int)vec.size(), name, "Data" );
			printf( "{" );
			for( size_t i = 0; i < vec.size(); i++ )
			{
				if( 0 == i % elementsPerLine )
					printf( "\n\t" );
				printElement( vec[ i ] );
				printf( ", " );
			}
			printf( "\n};" );
		}
	};

	template<bool signedValues>
	class Vec4Builder : public InnerIndexBase<uint32_t, Vector4i>
	{
		const char* storeTypeName() const override
		{
			return "uint32_t";
		}
		virtual uint32_t compressElement( const Vector4i& elt ) const override
		{
			uint32_t res = 0;
			for( int i = 3; i >= 0; i-- )
			{
				res = res << 8;
				int v = elt[ i ];
				if( signedValues )
				{
					if( v < -128 || v > +127 )
						__debugbreak();
					res |= (uint32_t)(uint8_t)v;
				}
				else
				{
					if( v < 0 || v > 0xFF )
						__debugbreak();
					res |= (uint32_t)v;
				}
			}
			return res;
		}
		virtual void printElement( uint32_t elt ) const override
		{
			printf( "0x%08X", (int)elt );
		}
	};

	class Vec2Builder : public InnerIndexBase<uint16_t, Vector2i>
	{
		const char* storeTypeName() const override
		{
			return "uint16_t";
		}
		virtual uint16_t compressElement( const Vector2i& elt ) const override
		{
			uint32_t res = 0;
			for( int i = 1; i >= 0; i-- )
			{
				int v = elt[ i ];
				if( v < 0 || v > 0xFF )
					__debugbreak();
				res = res << 8;
				res |= (uint32_t)v;
			}
			return res;
		}
		virtual void printElement( uint16_t elt ) const override
		{
			printf( "0x%04X", (int)elt );
		}
	};
	class Bool4Builder : public InnerIndexBase<uint8_t, std::array<bool, 4>>
	{
		const char* storeTypeName() const override
		{
			return "uint8_t";
		}
		virtual uint8_t compressElement( const std::array<bool, 4>& elt ) const override
		{
			uint8_t res = 0;
			if( elt[ 0 ] )
				res |= 1;
			if( elt[ 1 ] )
				res |= 2;
			if( elt[ 2 ] )
				res |= 4;
			if( elt[ 3 ] )
				res |= 8;
			return res;
		}
		virtual void printElement( uint8_t elt ) const override
		{
			uint32_t val = _pdep_u32( elt, 0x1111 );
			printf( "0b%04x", val );
		}
	};

	void printOuter( const char* name, const char* nameSuffix, const std::vector<uint32_t>& data )
	{
		constexpr size_t elementsPerLine = 32;
		const uint32_t maxValue = *data.rbegin();
		const char* dt;
		if( maxValue < 0x100 )
			dt = "uint8_t";
		else if( maxValue < 0x10000 )
			dt = "uint16_t";
		else
			dt = "uint32_t";

		printf( "static const std::array<%s, %i> %s%s = \n", dt, (int)data.size(), name, nameSuffix );
		printf( "{" );
		for( size_t i = 0; i < data.size(); i++ )
		{
			if( 0 == i % elementsPerLine )
				printf( "\n\t" );
			printf( "%i, ", (int)data[ i ] );
		}
		printf( "\n};\n" );
	}

	template<class E>
	using pfnGetVectorRef = const std::vector<std::vector<E>>& (*)( int i );

	template<class I, class E>
	void compressNestedVectors( int length, pfnGetVectorRef<E> pfn, const char* name )
	{
		I inner;
		std::vector<uint32_t> outer, outer2;
		for( int i = 0; i < length; i++ )
		{
			const std::vector<std::vector<E>>& vec2 = pfn( i );

			outer2.push_back( outer.size() );
			for( size_t j = 0; j < vec2.size(); j++ )
			{
				const std::vector<E>& vec1 = vec2[ j ];

				outer.push_back( inner.size() );
				for( const E& e : vec1 )
					inner.add( e );
			}
		}
		outer2.push_back( outer.size() );
		outer.push_back( inner.size() );

		printf( "\n" );
		printOuter( name, "Outer", outer2 );
		printOuter( name, "Inner", outer );
		inner.print( name );
	}
}  // namespace

void printCutTableData()
{
	g_file = _wfopen( LR"(C:\Temp\2remove\MeshRepair\CutTable2.inl)", L"w" );
	printf( "// This source file was generated automatically by a tool.\n// The tool is in C++, it consumed the data made by the code in auto_table.cpp" );
	compressNestedVectors<Vec4Builder<true>>( 64, &floatTetWild::CutTable::get_tet_confs, "TetraConfigs" );
	compressNestedVectors<Vec2Builder>( 64, &floatTetWild::CutTable::get_diag_confs, "DiagConfigs" );
	compressNestedVectors<Bool4Builder>( 64, &floatTetWild::CutTable::get_surface_conf, "SurfaceConfig" );
	compressNestedVectors<Vec4Builder<true>>( 64, &floatTetWild::CutTable::get_face_id_conf, "FaceIDs" );
	fclose( g_file );
}