#pragma once

namespace floatTetWild
{
	class VertexConnectedTets
	{
		std::vector<int> vec;
		friend void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result );
		friend void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result );

	  public:
		void add( int i )
		{
			vec.push_back( i );
		}
		void addRange( const std::vector<int>& other )
		{
			vec.insert( vec.end(), other.begin(), other.end() );
		}

		decltype( auto ) begin() const
		{
			return vec.begin();
		}
		decltype( auto ) end() const
		{
			return vec.end();
		}
		size_t size() const
		{
			return vec.size();
		}
		bool empty() const
		{
			return vec.empty();
		}
		bool remove( int idx )
		{
			auto it = std::find( vec.begin(), vec.end(), idx );
			if( it != vec.end() )
			{
				vec.erase( it );
				return true;
			}
			return false;
		}
		void eraseAt( size_t idx )
		{
			assert( idx < vec.size() );
			vec.erase( vec.begin() + idx );
		}
		void clear()
		{
			vec.clear();
		}
		int operator[]( size_t idx ) const
		{
			return vec[ idx ];
		}

		// Transform the elements v => map[ v ]
		void applyMapping( const std::vector<int>& map )
		{
			for( int& i : vec )
				i = map[ i ];
		}

		void sort()
		{
			std::sort( vec.begin(), vec.end() );
		}
	};

	void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result );
	void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result );
}  // namespace floatTetWild