#pragma once
#include "GeometricPrimitives.h"

namespace GEO2
{
	class NearestSearch
	{
	  public:
		void buildTree( const std::vector<vec3>& data );
		double getNearestPointSqDist( const vec3& pt ) const;
	};
}  // namespace GEO2