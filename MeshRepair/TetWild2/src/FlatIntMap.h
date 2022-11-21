#pragma once
#include <vector>

namespace floatTetWild
{
	// Implements API similar to std::map<int,int> on top of a sorted vector
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

		bool contains( int key ) const
		{
			auto it = std::lower_bound( vec.begin(), vec.end(), key, FindKey {} );
			if( it == vec.end() || it->key != key )
				return false;
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

		size_t size() const
		{
			return vec.size();
		}

		inline const std::pair<int, int>* begin() const
		{
			if( !vec.empty() )
				return (const std::pair<int, int>*)vec.data();
			return nullptr;
		}

		inline const std::pair<int, int>* end() const
		{
			if( !vec.empty() )
				return (const std::pair<int, int>*)( vec.data() + vec.size() );
			return nullptr;
		}
	};

	// Implements API similar to std::map<std::array<int,2>,int> on top of a sorted vector
	// In addition to that, the keys are treated as unsorted: [ 3, 7 ] and [ 7, 3 ] keys are equivalent,
	// i.e. this collection will resolve both keys into the same value
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

	// Despite the name, the sorting is optional.
	// Implemented for a weird use case in CutMesh class, where it first populates a sorted vector, then appends extra unsorted values to the end.
	class SortedIntSet
	{
		std::vector<int> vec;
#ifdef _DEBUG
		bool sorted = true;
#endif

	  public:
		void clear()
		{
			vec.clear();
#ifdef _DEBUG
			sorted = true;
#endif
		}

		void addSorted( int i )
		{
#ifdef _DEBUG
			assert( sorted );
#endif
			auto it = std::lower_bound( vec.begin(), vec.end(), i );
			if( it == vec.end() || *it != i )
				vec.insert( it, i );
		}
		void addUnsorted( int i )
		{
			vec.push_back( i );
#ifdef _DEBUG
			sorted = false;
#endif
		}

		size_t size() const
		{
			return vec.size();
		}
		int operator[]( size_t i ) const
		{
			return vec[ i ];
		}
	};
}  // namespace floatTetWild