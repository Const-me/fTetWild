#pragma once
#include "Geogram2.h"
#include "KdTree/kd_tree.h"

namespace GEO2
{
	class NearestSearch
	{
		BalancedKdTree tree;

	  public:

		NearestSearch();

		void buildTree( const std::vector<vec3>& data );

		double getNearestPointSqDist( const vec3& pt ) const;
	};
}  // namespace GEO2