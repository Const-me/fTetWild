#include <stdafx.h>
#include "TrackedSurfaceChanges.h"
using namespace floatTetWild;

void TrackedSurfaceChange::apply( std::array<std::vector<int>, 4>& rdi, const std::vector<std::array<std::vector<int>, 4>>& sourceVector ) const
{
	const std::array<std::vector<int>, 4>& rsi = sourceVector[ tet ];

	if( &rsi == &rdi )
	{
		// Applying in place
		std::array<std::vector<int>, 4> oldLanes;
		for( size_t i = 0; i < 4; i++ )
			oldLanes[ i ].swap( rdi[ i ] );

		for( size_t i = 0; i < 4; i++ )
		{
			const EntryChange& e = entries[ i ];

			if( e.sourceEntry < 4 )
			{
				rdi[ i ] = std::move( oldLanes[ e.sourceEntry ] );
				switch( e.flag )
				{
				case eNewItem::None:
					break;
				case eNewItem::Append:
					rdi[ i ].push_back( e.newItem );
					break;
				case eNewItem::Prepend:
					rdi[ i ].insert( rdi[ i ].begin(), e.newItem );
					break;
				}
				continue;
			}

			rdi[ i ].clear();
			if( e.flag != eNewItem::None )
				rdi[ i ].push_back( e.newItem );
		}
	}
	else
	{
		// Making a copy
		for( size_t i = 0; i < 4; i++ )
		{
			const EntryChange& e = entries[ i ];

			if( e.sourceEntry < 4 )
			{
				const std::vector<int>& source = rsi[ e.sourceEntry ];
				switch( e.flag )
				{
				case eNewItem::None:
					rdi[ i ] = source;
					break;
				case eNewItem::Append:
					rdi[ i ].reserve( source.size() + 1 );
					rdi[ i ] = source;
					rdi[ i ].push_back( e.newItem );
					break;
				case eNewItem::Prepend:
					rdi[ i ].clear();
					rdi[ i ].reserve( source.size() + 1 );
					rdi[ i ].push_back( e.newItem );
					rdi[ i ].insert( rdi[ i ].end(), source.begin(), source.end() );
					break;
				}
				continue;
			}

			rdi[ i ].clear();
			if( e.flag != eNewItem::None )
				rdi[ i ].push_back( e.newItem );
		}
	}
}