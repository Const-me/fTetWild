#include "GeogramDelaunay.h"
#include "geogram/delaunay/delaunay_3d.h"
#include "geogram/delaunay/parallel_delaunay_3d.h"

namespace
{
	template<class T>
	class Impl : public iDelaunay
	{
		T impl;

	  public:
		Impl()
			: impl( 3 )
		{
		}
		~Impl() override
		{
		}

		void compute( size_t count, const double* rsi ) override final
		{
			impl.set_vertices( (uint32_t)count, rsi );
		}

		size_t countElements() const override final
		{
			return impl.nb_cells();
		}

		const __m128i* getElements() const override final
		{
			return (const __m128i*)impl.cell_to_v();
		}
	};
}  // namespace

std::unique_ptr<iDelaunay> create( bool multithreaded )
{
	if( multithreaded )
		return std::make_unique<Impl<GEO::ParallelDelaunay3d>>();
	else
		return std::make_unique<Impl<GEO::Delaunay3d>>();
}