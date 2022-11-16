#pragma once
#include <atomic>
#include "logger.h"

// Accumulator to measure times of a specific function
class TimeMeasure
{
	std::atomic_size_t count = 0;
	std::atomic_uint64_t time = 0;

  public:
	// RAII class which measures it's own lifetime, i.e. time between constructor and destructor
	class MeasureRaii
	{
		TimeMeasure& source;
		uint64_t started;

	  public:
		MeasureRaii( TimeMeasure& tm );
		~MeasureRaii();
	};

	inline MeasureRaii measure()
	{
		return MeasureRaii { *this };
	}

	void logInfo( const Logger& log, double mulSeconds, const char* what ) const;
};

// A set of accumulators to measure some functions in this DLL
class TimeMeasures
{
	// Couple numbers to compute clock frequency of the measures
	const uint64_t constructedTsc, constructedQpc;

  public:
	TimeMeasures();

	// Print collected counters to the log
	void logInfo( const Logger& log ) const;

	// edge_collapsing_aux
	TimeMeasure edgeCollapsingAux; 
	// collapse_an_edge
	TimeMeasure collapseAnEdge;
	// insert_triangles_aux
	TimeMeasure insertTrianglesAux;
	// insert_one_triangle
	TimeMeasure insertOneTriangle;
	// find_cutting_tets
	TimeMeasure findCuttingTets;
	// edge_swapping
	TimeMeasure edgeSwapping;
	// subdivide_tets
	TimeMeasure subdivideTets;
};