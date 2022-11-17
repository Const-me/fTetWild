#include "stdafx.h"
#include "TrianglePartition.h"
#include "../Utils/AvxMath.h"
using namespace floatTetWild;

namespace
{
	// Utility class to load vertices from a vertex buffer, fast
	class VertexLoader
	{
		const Vector3* const rsi;
		const Vector3* const rsiLastAligned;
#ifdef _DEBUG
		const size_t vbLength;
#endif

	  public:
		VertexLoader( const std::vector<Vector3>& vb )
			: rsi( vb.data() )
			, rsiLastAligned( vb.data() + ( (ptrdiff_t)vb.size() - 1 ) )
#ifdef _DEBUG
			, vbLength( vb.size() )
#endif
		{
		}

		__forceinline __m256d load( int si ) const
		{
#ifdef _DEBUG
			assert( si >= 0 );
			assert( (uint32_t)si < vbLength );
#endif
			const Vector3* const p = rsi + (uint32_t)si;
			if( p < rsiLastAligned )
				return _mm256_loadu_pd( p->data() );
			else
				return AvxMath::loadDouble3( p->data() );
		}
	};

	static const double dblMax = DBL_MAX;

	inline void computeBox( const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d& i, __m256d& ax )
	{
		VertexLoader ldr { vb };

		i = _mm256_broadcast_sd( &dblMax );
		ax = AvxMath::vectorNegate( i );
		__m256d v;

		for( int fi : faces )
		{
			const Vector3i& indices = ib[ (uint32_t)fi ];

			v = ldr.load( indices[ 0 ] );
			i = _mm256_min_pd( i, v );
			ax = _mm256_max_pd( i, v );

			v = ldr.load( indices[ 1 ] );
			i = _mm256_min_pd( i, v );
			ax = _mm256_max_pd( i, v );

			v = ldr.load( indices[ 2 ] );
			i = _mm256_min_pd( i, v );
			ax = _mm256_max_pd( i, v );
		}
	}
}  // namespace

bool TrianglePartition::tryPartition(
  const std::vector<Vector3>& vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, __m256d clearance, std::array<std::vector<int>, 3>& result )
{
	if( faces.size() < 2 )
		return false;

	using namespace AvxMath;
	// Compute bounding box of all input faces
	__m256d boxMin, boxMax;
	computeBox( vb, ib, faces, boxMin, boxMax );

	// Find maximum dimension of that box, relative to the clearance
	const __m256d inputBoxSize = _mm256_sub_pd( boxMax, boxMin );
	const __m256d inputBoxRel = _mm256_div_pd( inputBoxSize, clearance );
	const __m256d max4 = vector3BroadcastMaximum( inputBoxRel );

	if( _mm256_cvtsd_f64( max4 ) < 2.0 )
		return false;  // Unlikely to succeed, the box is too small in relation to the clearance

	// Figure out direction of the partition
	const uint32_t splitCoordinateIndex = AvxMath::firstEqualLaneIndex( max4, inputBoxRel );
	assert( splitCoordinateIndex < 3 );

	// Call the scalar partitioning method, on the specified coordinate
	static_assert( sizeof( Vector3 ) == 3 * sizeof( double ) );
	const double* vbScalar = (const double*)( vb.data() ) + splitCoordinateIndex;
	const double c = vectorExtractLane( clearance, splitCoordinateIndex );
	return tryPartitionOnCoordinate( vbScalar, ib, faces, c, result );
}

bool TrianglePartition::tryPartitionOnCoordinate(
  const double* vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces, double clearance, std::array<std::vector<int>, 3>& result )
{
	buildEntries( vb, ib, faces );

	size_t leftEnd, middleEnd;
	if( !tryPartitionEntries( leftEnd, middleEnd, clearance ) )
		return false;

	makePartitions( leftEnd, middleEnd, result );
	return true;
}

inline void TrianglePartition::PartitionEntry::computeMinMax( const double* vb, const std::vector<Vector3i>& ib, int triIndex )
{
	const Vector3i& tri = ib[ triIndex ];

	size_t i0 = (size_t)(uint32_t)tri[ 0 ];
	size_t i1 = (size_t)(uint32_t)tri[ 1 ];
	size_t i2 = (size_t)(uint32_t)tri[ 2 ];

	i0 *= 3;
	i1 *= 3;
	i2 *= 3;

	const double c0 = vb[ i0 ];
	const double c1 = vb[ i1 ];
	const double c2 = vb[ i2 ];

	min = std::min( std::min( c0, c1 ), c2 );
	max = std::max( std::max( c0, c1 ), c2 );
}

struct TrianglePartition::SortEntries
{
	bool operator()( const PartitionEntry& a, const PartitionEntry& b )
	{
		return a.min < b.min;
	}
};

void TrianglePartition::buildEntries( const double* vb, const std::vector<Vector3i>& ib, const std::vector<int>& faces )
{
	// Initialize entries with min/max for every triangle in the input set
	entries.resize( faces.size() );
	for( size_t i = 0; i < faces.size(); i++ )
	{
		PartitionEntry& e = entries[ i ];
		e.computeMinMax( vb, ib, faces[ i ] );
		e.idxFaceInSet = (uint32_t)i;
		e.idxFaceInMesh = (uint32_t)faces[ i ];
	}

	// Sort entries my minimum coordinate
	std::sort( entries.begin(), entries.end(), SortEntries {} );

	// Initialize cumulative maximum of the entries
	double maxCumulative = entries[ 0 ].max;
	for( PartitionEntry& e : entries )
	{
		double ax = e.max;
		maxCumulative = std::max( maxCumulative, ax );
		e.max = maxCumulative;
	}
}

bool TrianglePartition::tryPartitionEntries( size_t& leftEnd, size_t& middleEnd, double clearance ) const
{
	double firstMax = entries[ 0 ].max;
	double lastMin = entries.back().min;
	if( firstMax + clearance >= lastMin )
	{
		// Clearance between first and last triangle is not enough
		return false;
	}

	// TODO: it's possible to accomplish the same with a binary search instead of linear one
	size_t mid1 = 0;
	size_t mid2 = entries.size() - 1;

	while( true )
	{
		mid1++;
		if( mid1 == mid2 )
			break;
		firstMax = entries[ mid1 ].max;
		if( firstMax + clearance >= lastMin )
			break;

		mid2--;
		if( mid1 == mid2 )
			break;
		lastMin = entries[ mid2 ].min;
		if( firstMax + clearance >= lastMin )
		{
			mid2++;
			break;
		}
	}
	leftEnd = mid1;
	middleEnd = mid2;
	return true;
}

void TrianglePartition::makePartitions( size_t leftEnd, size_t middleEnd, std::array<std::vector<int>, 3>& result )
{
	assert( leftEnd > 0 );
	assert( middleEnd >= leftEnd );
	assert( middleEnd < entries.size() );

	makePartition( 0, leftEnd, result[ 0 ] );
	makePartition( middleEnd, entries.size(), result[ 1 ] );
	if( middleEnd > leftEnd )
		makePartition( leftEnd, middleEnd, result[ 2 ] );
	else
	{
		// Unlikely, but this may happen - no triangles intersected the boundary between left and right partitions
		result[ 2 ].clear();
	}
}

struct TrianglePartition::SortTemp
{
	bool operator()( const EntryTemp& a, const EntryTemp& b )
	{
		return a.idxFaceInSet < b.idxFaceInSet;
	}
};

void TrianglePartition::makePartition( size_t begin, size_t end, std::vector<int>& dest )
{
	assert( end > begin );
	const size_t length = end - begin;

	if constexpr( preserveOrder )
	{
		entriesTemp.resize( length );
		for( size_t i = 0; i < length; i++ )
		{
			const PartitionEntry& rsi = entries[ i + begin ];
			EntryTemp& rdi = entriesTemp[ i ];

			// Copying two integers (8 bytes) with 1 instruction
			uint64_t val = *(const uint64_t*)( &rsi.idxFaceInSet );
			*(uint64_t*)( &rdi ) = val;
		}

		std::sort( entriesTemp.begin(), entriesTemp.end(), SortTemp {} );

		dest.resize( length );
		for( size_t i = 0; i < length; i++ )
			dest[ i ] = (int)entriesTemp[ i ].idxFaceInMesh;
	}
	else
	{
		dest.resize( length );
		for( size_t i = 0; i < length; i++ )
		{
			const PartitionEntry& rsi = entries[ i + begin ];
			dest[ i ] = (int)rsi.idxFaceInMesh;
		}
	}
}