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