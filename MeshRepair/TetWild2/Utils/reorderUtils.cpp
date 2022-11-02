#include "stdafx.h"
#include "reorderUtils.h"

void remapIndicesImpl( uint32_t* pointer, size_t count, const std::vector<uint32_t>& order )
{
	// Invert the order.
	// Geogram code does that in place, not sure I like their tradeoff
	// they're using sign bit to track things which limit length to 2G
	std::vector<uint32_t> inverted;
	const size_t length = order.size();
	inverted.resize( length );
	for( size_t i = 0; i < length; i++ )
		inverted[ order[ i ] ] = (uint32_t)i;

	// Remap input using the inverted ordering
	uint32_t* const pointerEnd = pointer + count;
	for( ; pointer < pointerEnd ; pointer++) 
	{
		uint32_t i = *pointer;
		i = inverted[ i ];
		*pointer = i;
	}
}