#pragma once

namespace floatTetWild
{
	// A specialized collection for the set of tetrahedral elements
	// Every vertex of the mesh has this collection, keeping indices of the tetra element connected to the vertex
	// The set is assumed to be unordered; various setIntersection() functions sort elements of the container, for performance reason
	class VertexConnectedTets
	{
		mutable std::vector<int> vec;
		mutable bool isSorted = false;
		friend void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result );
		friend void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result );

		void ensureSorted() const;

	  public:
		void add( int i )
		{
			vec.push_back( i );
			isSorted = false;
		}
		void addRange( const std::vector<int>& other )
		{
			vec.insert( vec.end(), other.begin(), other.end() );
			isSorted = false;
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
			if( isSorted )
			{
				// When sorted, we can use binary search to find elements
				// Note that when an element is found and removed, the vector remains sorted
				auto it = std::lower_bound( vec.begin(), vec.end(), idx );
				if( it != vec.end() && *it == idx )
				{
					vec.erase( it );
					return true;
				}
				return false;
			}
			else
			{
				auto it = std::find( vec.begin(), vec.end(), idx );
				if( it != vec.end() )
				{
					vec.erase( it );
					return true;
				}
				return false;
			}
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
			isSorted = false;
		}

		void sort()
		{
			ensureSorted();
		}
	};

	// Find elements present in both VertexConnectedTets collections, append to the output vector
	// The integers are added in the sorted order
	void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, std::vector<int>& result );

	// Find the elements present in all 3 VertexConnectedTets collections, append to the output vector
	// The integers are added in the sorted order
	void setIntersection( const VertexConnectedTets& a, const VertexConnectedTets& b, const VertexConnectedTets& c, std::vector<int>& result );

}  // namespace floatTetWild