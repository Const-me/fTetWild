#pragma once
#include <random>

namespace floatTetWild
{
	class Random
	{
	  public:
		static int next( int min, int max )
		{
			static std::mt19937 gen( 42 );
			int res = ( gen() - std::mt19937::min() ) / double( std::mt19937::max() - std::mt19937::min() ) * ( max - min ) + min;

			return res;
		}

		template<typename T>
		static void shuffle( std::vector<T>& vec )
		{
			for( int i = vec.size() - 1; i > 0; --i )
			{
				using std::swap;
				const int index = next( 0, i + 1 );
				swap( vec[ i ], vec[ index ] );
			}
		}
	};
}