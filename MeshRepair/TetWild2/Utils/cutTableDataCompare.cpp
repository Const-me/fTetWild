#include "stdafx.h"
#include "cutTableData.h"
#include "../src/auto_table.hpp"
#define CutTable zzzCutTable
#include "../src/CutTable2.h"
#undef CutTable

namespace
{
	template<class A, class B>
	void compareTableData( A a, B b )
	{
		for( int i = 0; i < 64; i++ )
		{
			const auto& i1 = a( i );
			// if( i == 63 ) __debugbreak();
			const auto& i2 = b( i );
			if( i1.size() != i2.size() )
				__debugbreak();

			for( size_t j = 0; j < i1.size(); j++ )
			{
				const auto& v1 = i1[ j ];
				const auto& v2 = i2[ j ];
				if( v1.size() != v2.size() )
					__debugbreak();
				for( size_t k = 0; k < v1.size(); k++ )
				{
					const auto e1 = v1[ k ];
					const auto e2 = v2[ k ];
					if( e1 != e2 )
						__debugbreak();
				}
			}
		}
	}
}  // namespace

void validateCutTableData()
{
	using namespace floatTetWild;

	compareTableData( []( int i ) { return CutTable::get_tet_confs( i ); }, []( int i ) { return CutTable2::get_tet_confs( i ); } );
	compareTableData( []( int i ) { return CutTable::get_diag_confs( i ); }, []( int i ) { return CutTable2::get_diag_confs( i ); } );
	compareTableData( []( int i ) { return CutTable::get_surface_conf( i ); }, []( int i ) { return CutTable2::get_surface_conf( i ); } );
	compareTableData( []( int i ) { return CutTable::get_face_id_conf( i ); }, []( int i ) { return CutTable2::get_face_id_conf( i ); } );
}