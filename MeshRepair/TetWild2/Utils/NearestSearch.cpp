#include "stdafx.h"
#include "NearestSearch.h"
using namespace GEO2;

NearestSearch::NearestSearch()
	: tree( 3 )
{
}

void NearestSearch::buildTree( const std::vector<vec3>& data )
{
	tree.set_points( (BalancedKdTree::index_t)data.size(), (const double*)data.data() );
}

double NearestSearch::getNearestPointSqDist( const vec3& pt ) const
{
	BalancedKdTree::index_t index;
	double res;
	tree.get_nearest_neighbors( 1, &pt.x, &index, &res );
	return res;
}