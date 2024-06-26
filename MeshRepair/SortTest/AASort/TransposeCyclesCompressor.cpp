﻿// You don’t need to compile this source file, unless you want to.
// See TransposePermutations.inl for the output data it generates.

// Compile with /FC (Full Path of Source Code File in Diagnostics)
// https://learn.microsoft.com/en-us/cpp/build/reference/fc-full-path-of-source-code-file-in-diagnostics?view=msvc-170
#include "AASort.h"
#if TRANSPOSE_CYCLE_COMPRESSOR

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <atlstr.h>
#include <atlpath.h>
#include <vector>
#include <stdio.h>
#include <assert.h>

namespace
{
	CString makeOutputPath()
	{
		CPath path;
		path.m_strPath = __FILE__;
		path.RemoveFileSpec();
		path.Append( L"TransposePermutations.inl" );
		return path.m_strPath;
	}

	template<class E>
	static void printVectorData( FILE* f, const std::vector<E>& src, size_t valuesPerLine )
	{
		for( size_t i = 0; i < src.size(); i++ )
		{
			if( 0 == ( i % valuesPerLine ) )
				fputs( "\n\t", f );
			fprintf( f, "%i, ", (int)src[ i ] );
		}
	}

	class CyclesCompressor
	{
		// Actual bytes there
		std::vector<uint8_t> inner;
		// 1 entry = one cycle, values are offsets into the inner vector
		std::vector<uint16_t> cycles;
		// 1 entry = complete permutation of the matrix, values are offsets into the cycles vector
		std::vector<uint16_t> permutations;

	  public:
		void addPermutation( const std::vector<std::vector<int>>& source )
		{
			permutations.push_back( cycles.size() );
			for( const auto& cycle : source )
			{
				cycles.push_back( inner.size() );
				for( int i : cycle )
				{
					assert( i >= 0 && i < 0x100 );
					inner.push_back( (uint8_t)i );
				}
			}
		}

		void finalize()
		{
			assert( cycles.size() < 0x10000 );
			permutations.push_back( (uint16_t)cycles.size() );

			assert( inner.size() < 0x10000 );
			cycles.push_back( (uint16_t)inner.size() );
		}

		void writeTheFile()
		{
			const CString path = makeOutputPath();
			FILE* f = nullptr;
			if( 0 != _wfopen_s( &f, path, L"w" ) || nullptr == f )
				throw std::exception( "fopen failed" );

			fprintf_s( f, "// This source file is generated by a tool. See TransposeCyclesCompressor.cpp for the source code\n" );

			fprintf_s( f, "static const std::array<uint8_t, %zu> s_permuteInner = \n{", inner.size() );
			printVectorData( f, inner, 32 );
			fputs( "\n};\n", f );

			fprintf_s( f, "static const std::array<uint16_t, %zu> s_permuteCycles = \n{", cycles.size() );
			printVectorData( f, cycles, 24 );
			fputs( "\n};\n", f );

			fprintf_s( f, "static const std::array<uint16_t, %zu> s_permuteOuter = \n{", permutations.size() );
			printVectorData( f, permutations, 24 );
			fputs( "\n};", f );

			fclose( f );
		}
	};

	class OuterCyclesTemp
	{
		// Length in vectors
		const size_t length;
		std::vector<bool> visited;
		std::vector<int> cycle;
		std::vector<std::vector<int>> allCycles;

		size_t elements() const
		{
			return length * 4;
		}

		size_t permute( size_t idx ) const
		{
			const size_t col = idx % 4;
			const size_t row = idx / 4;
			return col * length + row;
		}

		static void printCycle( const std::vector<int>& vec )
		{
			if( vec.size() <= 1 )
				return;
			printf( "%zu [ ", vec.size() );
			for( size_t j = 0; j < vec.size() - 1; j++ )
				printf( "%i, ", vec[ j ] );
			printf( "%i ]\n", *vec.rbegin() );
		}

	  public:
		OuterCyclesTemp( size_t len )
			: length( len )
		{
			visited.resize( len * 4, false );
		}

		void findAllCycles()
		{
			for( size_t i = 0; i < elements(); i++ )
			{
				if( visited[ i ] )
					continue;
				visited[ i ] = true;

				cycle.clear();
				cycle.push_back( (int)i );
				const size_t first = i;
				size_t next = permute( i );
				while( next != first )
				{
					visited[ next ] = true;
					cycle.push_back( (int)next );
					next = permute( next );
				}
				if( cycle.size() < 2 )
					continue;
				allCycles.push_back( cycle );
			}
		}

		void printAllCycles()
		{
			findAllCycles();
			for( const auto& c : allCycles )
				printCycle( c );
		}

		void compress( CyclesCompressor& comp )
		{
			comp.addPermutation( allCycles );
		}
	};

	void compressCycles( size_t sizeBlocks, CyclesCompressor& comp )
	{
		OuterCyclesTemp ct { sizeBlocks };
		printf( "=== Cycles for 4x%zu ===\n", sizeBlocks );
		// ct.printAllCycles();
		ct.findAllCycles();
		ct.compress( comp );
	}
}  // namespace

namespace AASort
{
	void compressTransposeCycles()
	{
		CyclesCompressor comp;
		const size_t count = 64;
		for( size_t i = 1; i <= count; i++ )
			compressCycles( i, comp );
		comp.finalize();
		comp.writeTheFile();
		__debugbreak();
	}
}  // namespace AASort

#endif	// TRANSPOSE_CYCLE_COMPRESSOR