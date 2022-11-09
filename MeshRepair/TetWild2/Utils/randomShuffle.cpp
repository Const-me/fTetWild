#include "stdafx.h"
#include "randomShuffle.h"
#include <random>

void randomShuffle( std::vector<int>& vec )
{
	// https://en.cppreference.com/w/cpp/algorithm/random_shuffle
	std::random_device rd;
	std::mt19937 g( rd() );
	std::shuffle( vec.begin(), vec.end(), g );
}