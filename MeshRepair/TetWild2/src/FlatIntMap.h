#pragma once
#include <vector>

namespace floatTetWild
{
	// Implements API similar to a subset of std::map<int,int> on top of a sorted vector
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

	// Implements API similar to a subset of std::map<std::array<int,2>,int> on top of a sorted vector
	// In addition to that, the keys are treated as unsorted: [ 3, 7 ] and [ 7, 3 ] keys are equivalent, they will resolve into the same value
	class FlatEdgeMap
	{
		struct alignas( 16 ) Entry
		{
			uint64_t key;
			int value;
		};
		std::vector<Entry> vec;

		struct FindKey
		{
			inline bool operator()( const Entry& e, uint64_t key ) const
			{
				return e.key < key;
			}
		};

		static inline uint64_t makeKey( int i1, int i2 )
		{
			uint64_t u1 = (uint32_t)i1;
			uint64_t u2 = (uint32_t)i2;
			uint64_t k = ( u1 << 32 ) | u2;
			uint64_t flipped = _rotr64( k, 32 );
			return std::max( k, flipped );
		}

	  public:
		void clear()
		{
			vec.clear();
		}

		void setAt( int i0, int i1, int value )
		{
			const uint64_t key = makeKey( i0, i1 );
			auto it = std::lower_bound( vec.begin(), vec.end(), key, FindKey {} );
			if( it == vec.end() || it->key != key )
				vec.insert( it, Entry { key, value } );
			else
				it->value = value;
		}

		void setAt( const std::array<int, 2>& arr, int value )
		{
			setAt( arr[ 0 ], arr[ 1 ], value );
		}

		bool tryLookup( int i0, int i1, int& val ) const
		{
			const uint64_t key = makeKey( i0, i1 );
			auto it = std::lower_bound( vec.begin(), vec.end(), key, FindKey {} );
			if( it == vec.end() || it->key != key )
				return false;
			val = it->value;
			return true;
		}

		bool tryLookup( const std::array<int, 2>& arr, int& val ) const
		{
			return tryLookup( arr[ 0 ], arr[ 1 ], val );
		}

		bool contains( int i0, int i1 ) const
		{
			const uint64_t key = makeKey( i0, i1 );
			auto it = std::lower_bound( vec.begin(), vec.end(), key, FindKey {} );
			if( it == vec.end() || it->key != key )
				return false;
			return true;
		}

		bool contains( const std::array<int, 2>& arr ) const
		{
			return contains( arr[ 0 ], arr[ 1 ] );
		}

		template<class Lambda>
		inline void enumerate( Lambda lambda ) const
		{
			for( const auto& e : vec )
			{
				const uint64_t& key = e.key;
				const std::array<int, 2>& arr = *(const std::array<int, 2>*)&key;
				lambda( arr );
			}
		}
	};
}  // namespace floatTetWild