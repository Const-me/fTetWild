#pragma once
#include <array>

namespace floatTetWild
{
	struct ElementInQueue
	{
		std::array<int, 2> v_ids;
		double weight;

		ElementInQueue() = default;
		ElementInQueue( const std::array<int, 2>& ids, double w )
			: v_ids( ids )
			, weight( w )
		{
		}
		ElementInQueue( int e0, int e1, double w )
			: v_ids { e0, e1 }
			, weight( w )
		{
		}
	};

	struct cmp_l
	{
		bool operator()( const ElementInQueue& e1, const ElementInQueue& e2 )
		{
			if( e1.weight == e2.weight )
				return e1.v_ids > e2.v_ids;
			return e1.weight < e2.weight;
		}
	};

	struct cmp_s
	{
		bool operator()( const ElementInQueue& e1, const ElementInQueue& e2 )
		{
			if( e1.weight == e2.weight )
				return e1.v_ids < e2.v_ids;
			return e1.weight > e2.weight;
		}
	};
}