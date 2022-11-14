#include <iostream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include "AASort/AASort.h"

std::vector<int> createTest()
{
	srand( 0 );
	constexpr int len = 9999;
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
	std::vector<int> sortedStd, sortedAA;
	sortStd( src, sortedStd );
	AASort::sortVector( src, sortedAA );
	if( sortedStd != sortedAA )
		__debugbreak();
	std::vector<int> sortedAAInplace = src;
	AASort::sortVector( sortedAAInplace );
	if( sortedStd != sortedAAInplace )
		__debugbreak();

	std::cout << "Hello World!\n";
}