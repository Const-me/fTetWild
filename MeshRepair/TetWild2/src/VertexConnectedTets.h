#pragma once

#define CONN_TETS_SORTED_COPY 0

namespace floatTetWild
{
	class VertexConnectedTets
	{
		std::vector<int> vec;
		friend void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result );
		friend void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result );
#if CONN_TETS_SORTED_COPY
		mutable std::vector<int> sortedVec;
		const std::vector<int>& makeSortedVector() const;
		void flushSortedCopy()
		{
			sortedVec.clear();
		}
#else
		inline void flushSortedCopy()
		{
		}
#endif

	  public:
		void add( int i )
		{
			vec.push_back( i );
			flushSortedCopy();
		}
		void addRange( const std::vector<int>& other )
		{
			vec.insert( vec.end(), other.begin(), other.end() );
			flushSortedCopy();
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
				flushSortedCopy();
				return true;
			}
			return false;
		}
		void eraseAt( size_t idx )
		{
			assert( idx < vec.size() );
			vec.erase( vec.begin() + idx );
			flushSortedCopy();
		}
		void clear()
		{
			vec.clear();
			flushSortedCopy();
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
			flushSortedCopy();
		}

		void sort()
		{
			std::sort( vec.begin(), vec.end() );
#if CONN_TETS_SORTED_COPY
			sortedVec = vec;
#endif
		}
	};

	void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result );
	void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result );
}  // namespace floatTetWild