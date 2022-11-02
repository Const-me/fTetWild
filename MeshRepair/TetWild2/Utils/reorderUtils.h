#pragma once

template<class E, class TA>
inline void reorderVector( std::vector<E, TA>& vec, const std::vector<uint32_t>& order )
{
	const size_t length = vec.size();
	assert( length == order.size() );

	std::vector<E, TA> newVec;
	newVec.resize( length );

	for( size_t i = 0; i < length; i++ )
		newVec[ i ] = vec[ order[ i ] ];
	newVec.swap( vec );
}

void remapIndicesImpl( uint32_t* pointer, size_t length, const std::vector<uint32_t>& order );

template<class E, class TA>
inline void remapIndices( std::vector<E, TA>& vec, const std::vector<uint32_t>& order )
{
	assert( 0 == ( sizeof( E ) % 4 ) );

	size_t length = vec.size() * ( sizeof( E ) / 4 );
	remapIndicesImpl( (uint32_t*)vec.data(), length, order );
}