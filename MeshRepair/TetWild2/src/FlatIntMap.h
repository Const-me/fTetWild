#pragma once
#include <vector>

namespace floatTetWild
{
	// This class implements API similar to a subset of std::map<int,int> on top of a sorted vector
	// The elements are at sequential addresses, this structure is cache friendly, and when reused it retains the memory.
	class FlatIntMap
	{
		struct alignas( 8 ) Entry
		{
			int key, value;
		};
		std::vector<Entry> vec;

		struct FindKey
		{
			inline bool operator()( const Entry& e, int key ) const
			{
				return e.key < key;
			}
		};

	  public:
		void clear()
		{
			vec.clear();
		}

		void setAt( int key, int value )
		{
			auto it = std::lower_bound( vec.begin(), vec.end(), key, FindKey {} );
			if( it == vec.end() || it->key != key )
				vec.insert( it, Entry { key, value } );
			else
				it->value = value;
		}

		bool tryLookup( int key, int& val ) const
		{
			auto it = std::lower_bound( vec.begin(), vec.end(), key, FindKey {} );
			if( it == vec.end() || it->key != key )
				return false;
			val = it->value;
			return true;
		}

		int operator[]( int key ) const
		{
			auto it = std::lower_bound( vec.begin(), vec.end(), key, FindKey {} );
			if( it != vec.end() && it->key == key )
				return it->value;
			throw std::exception( "Key not found" );
		}

		// Inserts a new element into the container, if there is no element with the key in the container.
		// Returns a pair consisting of a pointer to the inserted value, or the already-existing value if no insertion happened, and a bool denoting
		// whether the insertion took place (true if insertion happened, false if it did not).
		std::pair<int*, bool> emplace( int key )
		{
			auto it = std::lower_bound( vec.begin(), vec.end(), key, FindKey {} );
			if( it != vec.end() && it->key == key )
				return std::make_pair( &it->value, false );

			it = vec.insert( it, Entry { key, -1 } );
			return std::make_pair( &it->value, true );
		}
	};
}  // namespace floatTetWild