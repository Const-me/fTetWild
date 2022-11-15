#pragma once
#include <array>

namespace floatTetWild
{
	// Small buffer of elements on stack. It may contain any count of element within [ 0 .. capacity ]
	template<class E, size_t capacity>
	class SmallBuffer
	{
		size_t length = 0;
		std::array<E, capacity> arr;

	  public:
		size_t size() const
		{
			return length;
		}
		bool empty() const
		{
			return length == 0;
		}
		void push_back( const E& element )
		{
			arr[ length ] = element;
			length++;
		}
		E& emplace_back()
		{
			E& ref = arr[ length ];
			length++;
			return ref;
		}

		const E& operator[]( size_t index ) const
		{
			return arr[ index ];
		}

		decltype( auto ) begin()
		{
			return arr.begin();
		}

		decltype( auto ) end()
		{
			return arr.begin() + length;
		}

		const E* data() const
		{
			assert( !empty() );
			return arr.data();
		}
	};
}  // namespace floatTetWild