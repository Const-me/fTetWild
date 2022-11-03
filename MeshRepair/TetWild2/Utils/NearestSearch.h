#pragma once
#include "GeometricPrimitives.h"

namespace GEO2
{
	class NearestSearch
	{
	  public:
		void buildTree( size_t count, const vec3* data );
		double getNearestPointSqDist( const vec3& pt ) const;
	};
}  // namespace GEO2