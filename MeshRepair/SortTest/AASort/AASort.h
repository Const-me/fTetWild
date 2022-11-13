#pragma once
#include <stdint.h>
#include <vector>

// AA-Sort: A New Parallel Sorting Algorithm for Multi-Core SIMD Processors
// Hiroshi Inoue, Takao Moriyama, Hideaki Komatsu and Toshio Nakatani, IBM Tokyo Research Laboratory
// DOI: 10.1109/PACT.2007.4336211
// This particular implementation requires at least SSE 4.1, AVX is recommended
namespace AASort
{
	// Sort the vector of signed 32-bit integers
	void sortVector( const std::vector<int>& source, std::vector<int>& dest );

	// Sort the vector of unsigned 32-bit integers
	void sortVector( const std::vector<uint32_t>& source, std::vector<uint32_t>& dest );
};