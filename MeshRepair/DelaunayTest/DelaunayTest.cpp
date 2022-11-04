#include "pch.h"
#include <iostream>
#include <vector>
#include <stdlib.h>
#include "../GeogramDelaunay/GeogramDelaunay.h"

int main()
{
	std::unique_ptr<iDelaunay> up = iDelaunay::create( false );

	std::vector<double> arr;
	const size_t countPoints = 1000;
	for( size_t i = 0; i < countPoints * 3; i++ )
		arr.push_back( rand() );
	std::cout << "Starting..\n";
	up->compute( countPoints, arr.data() );
	std::cout << "Completed, " << up->countElements() << " tetrahedra\n";
}