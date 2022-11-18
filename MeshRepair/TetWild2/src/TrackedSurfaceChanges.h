#pragma once
#include <vector>
#include <array>

namespace floatTetWild
{
	// A single change to the tracked surface collection
	struct TrackedSurfaceChange
	{
		int tet;

		enum struct eNewItem : uint8_t
		{
			None = 0,
			Append = 1,
			Prepend = 2
		};

		struct EntryChange
		{
			// When in [ 0 .. 3 ] interval, source index of the corresponding element.
			// When 0xFF, the vector is empty
			uint8_t sourceEntry;
			eNewItem flag;
			int newItem;
		};

		std::array<EntryChange, 4> entries;

		void makeCopyWithAppend( int id, uint8_t lane, int val )
		{
			tet = id;
			for( size_t i = 0; i < 4; i++ )
			{
				entries[ i ].sourceEntry = (uint8_t)i;
				entries[ i ].flag = eNewItem::None;
				entries[ i ].newItem = -1;
			}

			entries[ lane ].flag = eNewItem::Append;
			entries[ lane ].newItem = val;
		}

		void makeEmpty( int id )
		{
			tet = id;
			for( auto& e : entries )
			{
				e.sourceEntry = 0xFF;
				e.flag = eNewItem::None;
				e.newItem = -1;
			}
		}

		void prepend( uint8_t lane, int item )
		{
			auto& e = entries[ lane ];
			e.flag = eNewItem::Prepend;
			e.newItem = item;
		}

		void setSourceCopy( uint8_t lane, uint8_t sourceLane )
		{
			entries[ lane ].sourceEntry = sourceLane;
		}

		void apply( std::array<std::vector<int>, 4>& rdi, const std::vector<std::array<std::vector<int>, 4>>& rsi ) const;
	};

	// A set of changes to the tracked surface collection
	using TSChanges = std::vector<TrackedSurfaceChange>;
}  // namespace floatTetWild