#pragma once
#include <memory>
#include <emmintrin.h>

class iDelaunay
{
  public:
	virtual ~iDelaunay()
	{
	}

	virtual void compute( size_t count, const double* rsi ) = 0;

	virtual size_t countElements() const = 0;

	virtual const __m128i* getElements() const = 0;

	static std::unique_ptr<iDelaunay> create( bool multithreaded );
};