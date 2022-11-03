#include "stdafx.h"
#include "TetrahedralMeshAABB.h"
#include "TetrahedralMesh.h"
#include "../BoundingBox.hpp"
#include <geogram/numerics/predicates.h>

namespace
{
	/**
	 * \brief Computes the maximum node index in a subtree
	 * \param[in] node_index node index of the root of the subtree
	 * \param[in] b first facet index in the subtree
	 * \param[in] e one position past the last facet index in the subtree
	 * \return the maximum node index in the subtree rooted at \p node_index
	 */
	size_t max_node_index( size_t node_index, size_t b, size_t e )
	{
		assert( e > b );
		if( b + 1 == e )
			return node_index;

		size_t m = b + ( e - b ) / 2;
		size_t childl = 2 * node_index;
		size_t childr = 2 * node_index + 1;
		return std::max( max_node_index( childl, b, m ), max_node_index( childr, m, e ) );
	}
}  // namespace

TetrahedralMeshAABB::TetrahedralMeshAABB( TetrahedralMesh& sourceMesh, bool reorder )
	: mesh( sourceMesh )
{
	if( reorder )
		sourceMesh.reorderMorton();

	const size_t nb = sourceMesh.countElements();
	// +1 because size == max_index + 1
	boxes.resize( max_node_index( 1, 0, nb ) + 1 );
	initializeRecursive( 1, 0, nb );
}

namespace
{
	inline const double* vertexPointer( const TetrahedralMesh& mesh, int i )
	{
		return &mesh.getVertex( i ).x;
	}
}  // namespace

void TetrahedralMeshAABB::computeBox( BBox& rdi, size_t index ) const
{
	const __m128i e = mesh.getElement( index );

	BoundingBox box { vertexPointer( mesh, _mm_cvtsi128_si32( e ) ) };
	box.extend( vertexPointer( mesh, _mm_extract_epi32( e, 1 ) ) );
	box.extend( vertexPointer( mesh, _mm_extract_epi32( e, 2 ) ) );
	box.extend( vertexPointer( mesh, _mm_extract_epi32( e, 3 ) ) );

	box.store( rdi.min.data(), rdi.max.data() );
}

inline void TetrahedralMeshAABB::mergeBoxes( BBox& rdi, const BBox& a, const BBox& b )
{
	// min.xy
	__m128d v = _mm_load_pd( a.min.data() );
	v = _mm_min_pd( v, _mm_load_pd( b.min.data() ) );
	_mm_store_pd( rdi.min.data(), v );
	// min.z
	rdi.min[ 2 ] = std::min( a.min[ 2 ], b.min[ 2 ] );
	// max.x
	rdi.max[ 0 ] = std::max( a.max[ 0 ], b.max[ 0 ] );
	// max.yz - note it's aligned by 16 bytes just like min.xy
	v = _mm_load_pd( &a.max[ 1 ] );
	v = _mm_max_pd( v, _mm_load_pd( &b.max[ 1 ] ) );
	_mm_store_pd( &rdi.max[ 1 ], v );
}

void TetrahedralMeshAABB::initializeRecursive( size_t node_index, size_t b, size_t e )
{
	assert( node_index < boxes.size() );
	assert( b != e );
	if( b + 1 == e )
	{
		computeBox( boxes[ node_index ], b );
		return;
	}

	size_t m = b + ( e - b ) / 2;
	size_t childl = 2 * node_index;
	size_t childr = 2 * node_index + 1;
	assert( childl < boxes.size() );
	assert( childr < boxes.size() );
	initializeRecursive( childl, b, m );
	initializeRecursive( childr, m, e );

	assert( childl < boxes.size() );
	assert( childr < boxes.size() );
	mergeBoxes( boxes[ node_index ], boxes[ childl ], boxes[ childr ] );
}

inline bool TetrahedralMeshAABB::BBox::containsPoint( const GEO::vec3& p ) const
{
	const __m128d p1 = _mm_loadu_pd( &p.x );
	const __m128d p2 = _mm_load_sd( &p.z );

	__m128d v = _mm_load_pd( &min[ 0 ] );
	__m128d oob = _mm_cmplt_pd( p1, v );  // p.xy < min.xy
	v = _mm_load_sd( &min[ 2 ] );
	oob = _mm_or_pd( oob, _mm_cmplt_pd( p2, v ) );	// |= p.z < min.z

	v = _mm_loadu_pd( &max[ 0 ] );
	oob = _mm_or_pd( oob, _mm_cmpgt_pd( p1, v ) );	// |= p.xy > max.xy
	v = _mm_load_sd( &max[ 2 ] );
	oob = _mm_or_pd( oob, _mm_cmpgt_pd( p2, v ) );	// |= p.z > max.z

	return (bool)_mm_testz_pd( oob, oob );
}

namespace
{
	inline bool elementContainsPoint( const TetrahedralMesh& mesh, size_t elt, const GEO::vec3& p )
	{
		const __m128i indices = mesh.getElement( elt );

		using namespace GEO;
		const vec3& p0 = mesh.getVertex( (uint32_t)_mm_cvtsi128_si32( indices ) );
		const vec3& p1 = mesh.getVertex( (uint32_t)_mm_extract_epi32( indices, 1 ) );
		const vec3& p2 = mesh.getVertex( (uint32_t)_mm_extract_epi32( indices, 2 ) );
		const vec3& p3 = mesh.getVertex( (uint32_t)_mm_extract_epi32( indices, 3 ) );

		Sign s[ 4 ];
		s[ 0 ] = PCK::orient_3d( p, p1, p2, p3 );
		s[ 1 ] = PCK::orient_3d( p0, p, p2, p3 );
		s[ 2 ] = PCK::orient_3d( p0, p1, p, p3 );
		s[ 3 ] = PCK::orient_3d( p0, p1, p2, p );

		return ( ( s[ 0 ] >= 0 && s[ 1 ] >= 0 && s[ 2 ] >= 0 && s[ 3 ] >= 0 ) || ( s[ 0 ] <= 0 && s[ 1 ] <= 0 && s[ 2 ] <= 0 && s[ 3 ] <= 0 ) );
	}
}  // namespace

uint32_t TetrahedralMeshAABB::containingElement( const GEO::vec3& p, size_t n, size_t b, size_t e ) const
{
	if( !boxes[ n ].containsPoint( p ) )
		return NO_TET;

	if( e == b + 1 )
	{
		if( elementContainsPoint( mesh, b, p ) )
			return b;
		else
			return NO_TET;
	}

	size_t m = b + ( e - b ) / 2;
	size_t childl = 2 * n;
	size_t childr = 2 * n + 1;

	uint32_t result = containingElement( p, childl, b, m );
	if( result == NO_TET )
		result = containingElement( p, childr, m, e );
	return result;
}

uint32_t TetrahedralMeshAABB::containingElement( const GEO::vec3& p ) const
{
	return containingElement( p, 1, 0, mesh.countElements() );
}