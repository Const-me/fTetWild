#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include "AASort/AASort.h"
#include "AASort/OuterCyclesCompressor.h"

std::vector<int> createTest()
{
	srand( 0 );
	constexpr int len = 1777;
	std::vector<int> res;
	res.resize( len );
	for( int& i : res )
		i = rand();
	return res;
}

void sortStd( const std::vector<int>& source, std::vector<int>& dest )
{
	dest = source;
	std::sort( dest.begin(), dest.end() );
}

int main()
{
	// AASort::compressTransposeCycles();

	std::vector<int> src = createTest();
	std::vector<int> sortedStd, sorgedAA;
	sortStd( src, sortedStd );
	AASort::sortVector( src, sorgedAA );
	if( sortedStd != sorgedAA )
		__debugbreak();

	std::cout << "Hello World!\n";
}