#pragma once

class EdgesSet
{
	std::vector<std::array<int, 2>> edges;

  public:
	void clear()
	{
		edges.clear();
	}

	template<class Lambda>
	inline void enumerate( Lambda lambda ) const
	{
		for( const auto& e : edges )
			lambda( e[ 0 ], e[ 1 ] );
	}

	inline void add( int e0, int e1 )
	{
		edges.push_back( { { e0, e1 } } );
	}

	inline void addSorted( int e0, int e1 )
	{
		if( e0 > e1 )
			std::swap( e0, e1 );
		add( e0, e1 );
	}

	inline void reserve( size_t len )
	{
		edges.reserve( len );
	}

	void sortUnique();

	size_t size() const
	{
		return edges.size();
	}
};

// Sorted collection of integer edge IDs.
// Bad fit for keeping thousands of elements (inserting becomes expensive), good fit for local operations with small count of elements
class SortedEdgesSet
{
	std::vector<uint64_t> vec;

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

	void add( int i0, int i1 )
	{
		const uint64_t key = makeKey( i0, i1 );
		auto it = std::lower_bound( vec.begin(), vec.end(), key );
		if( it == vec.end() || *it != key )
			vec.insert( it, key );
	}

	size_t size() const
	{
		return vec.size();
	}

	template<class Lambda>
	void enumerate( Lambda lambda ) const
	{
		for( const uint64_t& u64 : vec )
		{
			const std::array<int, 2>& arr = *(const std::array<int, 2>*)&u64;
			lambda( arr );
		}
	}

	const std::array<int, 2>& operator[]( size_t index ) const
	{
		const uint64_t& u64 = vec[ index ];
		return *(const std::array<int, 2>*)&u64;
	}
};