#include "stdafx.h"
#include "LocalOperations.h"
#include "LocalOperations2.h"
#include <Utils/lowLevel.h>

namespace
{
	inline double pow2( double x )
	{
		return x * x;
	}

	inline double cubicRoot( double x )
	{
		return floatTetWild::cbrt( x );
	}

	constexpr double magic1 = 0.577350269189626;  // sqrt(1/3)
	constexpr double magic2 = 1.15470053837925;	  // sqrt(4/3)
	constexpr double magic3 = 0.408248290463863;  // sqrt(1/6)
	constexpr double magic4 = 1.22474487139159;	  // sqrt(3/2)
	constexpr double magic5 = 0.707106781186548;  // sqrt(0.5)
}

void floatTetWild::AMIPS_hessian_v2( const std::array<double, 12>& arr, Matrix3& result_0 )
{
	const double v0_x = arr[ 0 ];
	const double v0_y = arr[ 1 ];
	const double v0_z = arr[ 2 ];
	const double v1_x = arr[ 3 ];
	const double v1_y = arr[ 4 ];
	const double v1_z = arr[ 5 ];
	const double v2_x = arr[ 6 ];
	const double v2_y = arr[ 7 ];
	const double v2_z = arr[ 8 ];
	const double v3_x = arr[ 9 ];
	const double v3_y = arr[ 10 ];
	const double v3_z = arr[ 11 ];

	const double helper_5 = magic1 * v0_x;
	const double helper_21 = magic1 * v0_y;
	const double helper_42 = magic1 * v0_z;

	const double helper_9 = magic1 * v3_x;
	const double helper_23 = magic1 * v3_y;
	const double helper_44 = magic1 * v3_z;

	const double helper_7 = magic2 * v1_x;
	const double helper_22 = magic2 * v1_y;
	const double helper_43 = magic2 * v1_z;

	const double helper_12 = magic3 * v0_y;
	const double helper_14 = magic3 * v1_y;
	const double helper_18 = magic3 * v3_y;
	const double helper_25 = magic3 * v0_x;
	const double helper_26 = magic3 * v1_x;
	const double helper_29 = magic3 * v3_x;
	const double helper_34 = magic3 * v0_z;
	const double helper_36 = magic3 * v1_z;
	const double helper_39 = magic3 * v3_z;

	const double helper_28 = magic4 * v2_x;
	const double helper_16 = magic4 * v2_y;
	const double helper_38 = magic4 * v2_z;

	const double helper_3 = v0_z - v3_z;
	const double helper_10 = helper_5 - helper_7 + helper_9;
	const double helper_19 = helper_12 + helper_14 - helper_16 + helper_18;
	const double helper_20 = helper_10 * helper_19;
	const double helper_24 = helper_21 - helper_22 + helper_23;
	const double helper_30 = helper_25 + helper_26 - helper_28 + helper_29;
	const double helper_31 = helper_24 * helper_30;
	const double helper_32 = helper_3 * ( helper_20 - helper_31 );
	const double helper_33 = v0_x - v3_x;
	const double helper_40 = helper_34 + helper_36 - helper_38 + helper_39;
	const double helper_41 = helper_24 * helper_40;
	const double helper_45 = helper_42 - helper_43 + helper_44;
	const double helper_46 = helper_19 * helper_45;
	const double helper_47 = helper_41 - helper_46;
	const double helper_48 = helper_33 * helper_47;
	const double helper_49 = v0_y - v3_y;
	const double helper_50 = helper_10 * helper_40;
	const double helper_51 = helper_30 * helper_45;
	const double helper_52 = helper_50 - helper_51;
	const double helper_53 = helper_49 * helper_52;
	const double helper_54 = helper_32 + helper_48 - helper_53;
	const double helper_55 = pow2( helper_54 );
	const double helper_56 = 1.0 / cubicRoot( helper_55 );
	const double helper_57 = 1.0 * v2_x - 3.0 * v0_x + 1.0 * v1_x + 1.0 * v3_x;
	const double helper_58 = magic5 * v1_y;
	const double helper_59 = magic5 * v2_y;
	const double helper_60 = helper_58 - helper_59;
	const double helper_61 = helper_3 * helper_60;
	const double helper_62 = magic5 * v1_z - magic5 * v2_z;
	const double helper_63 = helper_49 * helper_62;
	const double helper_64 = helper_47 + helper_61 - helper_63;
	const double helper_65 = 1.33333333333333 / helper_54;
	const double helper_66 = 1.0 / helper_55;
	const double helper_67 = 0.5 * v2_x + 0.5 * v1_x;
	const double helper_68 = -1.5 * v0_x + helper_67 + 0.5 * v3_x;
	const double helper_69 = 0.5 * v0_x + helper_67 - 1.5 * v3_x;
	const double helper_70 = -1.5 * v2_x + 0.5 * v0_x + 0.5 * v1_x + 0.5 * v3_x;
	const double helper_71 = 0.5 * v2_x + 0.5 * v0_x - 1.5 * v1_x + 0.5 * v3_x;
	const double helper_72 = 0.5 * v1_y + 0.5 * v2_y;
	const double helper_73 = -1.5 * v0_y + 0.5 * v3_y + helper_72;
	const double helper_74 = 0.5 * v0_y - 1.5 * v3_y + helper_72;
	const double helper_75 = 0.5 * v0_y + 0.5 * v1_y - 1.5 * v2_y + 0.5 * v3_y;
	const double helper_76 = 0.5 * v0_y - 1.5 * v1_y + 0.5 * v2_y + 0.5 * v3_y;
	const double helper_77 = 0.5 * v1_z + 0.5 * v2_z;
	const double helper_78 = -1.5 * v0_z + 0.5 * v3_z + helper_77;
	const double helper_79 = 0.5 * v0_z - 1.5 * v3_z + helper_77;
	const double helper_80 = 0.5 * v0_z + 0.5 * v3_z + 0.5 * v1_z - 1.5 * v2_z;
	const double helper_81 = 0.5 * v0_z + 0.5 * v3_z - 1.5 * v1_z + 0.5 * v2_z;
	const double helper_82 = v0_z * helper_78 + v0_y * helper_73 + v1_y * helper_76 + v2_y * helper_75 + v3_y * helper_74 +
							 v3_z * helper_79 + v2_x * helper_70 + v1_z * helper_81 + v2_z * helper_80 + v0_x * helper_68 +
							 v1_x * helper_71 + helper_69 * v3_x;
	const double helper_83 = 0.444444444444444 * helper_66 * helper_82;
	const double helper_84 = helper_66 * helper_82;
	const double helper_85 = -helper_32 - helper_48 + helper_53;
	const double helper_86 = 1.0 / helper_85;
	const double helper_87 = helper_86 / cubicRoot( pow2( helper_85 ) );
	const double helper_88 = magic5 * v1_x;
	const double helper_89 = magic5 * v2_x;
	const double helper_90 = helper_88 - helper_89;
	const double helper_91 = 0.666666666666667 * helper_10 * helper_40 + 0.666666666666667 * helper_3 * helper_90 - 0.666666666666667 * helper_30 * helper_45 -
							 0.666666666666667 * helper_33 * helper_62;
	const double helper_92 = -3.0 * v0_y + 1.0 * v1_y + 1.0 * v2_y + 1.0 * v3_y;
	const double helper_93 = -v0_y + v3_y;
	const double helper_94 = -v0_z + v3_z;
	const double helper_95 = -helper_21 + helper_22 - helper_23;
	const double helper_96 = -helper_34 - helper_36 + helper_38 - helper_39;
	const double helper_97 = -helper_42 + helper_43 - helper_44;
	const double helper_98 = -helper_12 - helper_14 + helper_16 - helper_18;
	const double helper_99 = -0.666666666666667 * helper_60 * helper_94 + 0.666666666666667 * helper_62 * helper_93 +
							 0.666666666666667 * helper_95 * helper_96 - 0.666666666666667 * helper_97 * helper_98;
	const double helper_100 = helper_3 * helper_90;
	const double helper_101 = helper_33 * helper_62;
	const double helper_102 = helper_100 - helper_101 + helper_52;
	const double helper_103 = -helper_60 * helper_94 + helper_62 * helper_93 + helper_95 * helper_96 - helper_97 * helper_98;
	const double helper_104 = 0.444444444444444 * helper_102 * helper_103 * helper_82 * helper_86 + helper_57 * helper_91 - helper_92 * helper_99;
	const double helper_105 =
	  1.85037170770859e-17 * v0_z * helper_78 + 1.85037170770859e-17 * v0_y * helper_73 + 1.85037170770859e-17 * v1_y * helper_76 +
	  1.85037170770859e-17 * v2_y * helper_75 + 1.85037170770859e-17 * v3_y * helper_74 + 1.85037170770859e-17 * v3_z * helper_79 +
	  1.85037170770859e-17 * v2_x * helper_70 + 1.85037170770859e-17 * v1_z * helper_81 + 1.85037170770859e-17 * v2_z * helper_80 +
	  1.85037170770859e-17 * v0_x * helper_68 + 1.85037170770859e-17 * v1_x * helper_71 + 1.85037170770859e-17 * helper_69 * v3_x;
	const double helper_106 = helper_64 * helper_82 * helper_86;
	const double helper_107 = -0.666666666666667 * helper_10 * helper_19 + 0.666666666666667 * helper_24 * helper_30 +
							  0.666666666666667 * helper_33 * helper_60 - 0.666666666666667 * helper_49 * helper_90;
	const double helper_108 = -3.0 * v0_z + 1.0 * v3_z + 1.0 * v1_z + 1.0 * v2_z;
	const double helper_109 = -helper_20 + helper_31 + helper_33 * helper_60 - helper_49 * helper_90;
	const double helper_110 = 0.444444444444444 * helper_109 * helper_82 * helper_86;
	const double helper_111 = helper_103 * helper_110 + helper_107 * helper_57 - helper_108 * helper_99;
	const double helper_112 = -v0_x + v3_x;
	const double helper_113 = -helper_88 + helper_89;
	const double helper_114 = -helper_5 + helper_7 - helper_9;
	const double helper_115 = -helper_25 - helper_26 + helper_28 - helper_29;
	const double helper_116 = helper_82 * helper_86 * ( helper_112 * helper_62 + helper_113 * helper_94 + helper_114 * helper_96 - helper_115 * helper_97 );
	const double helper_117 = -helper_100 + helper_101 - helper_50 + helper_51;
	const double helper_118 = -helper_102 * helper_110 + helper_107 * helper_92 + helper_108 * helper_91;
	const double helper_119 =
	  helper_82 * helper_86 * ( helper_112 * ( -helper_58 + helper_59 ) - helper_113 * helper_93 - helper_114 * helper_98 + helper_115 * helper_95 );
	result_0( 0, 0 ) = helper_56 * ( helper_57 * helper_64 * helper_65 - pow2( helper_64 ) * helper_83 +
									 0.666666666666667 * helper_64 * helper_84 * ( -helper_41 + helper_46 - helper_61 + helper_63 ) + 3.0 );
	result_0( 0, 1 ) = helper_87 * ( helper_104 - helper_105 * v1_z + helper_106 * helper_91 );
	result_0( 0, 2 ) = helper_87 * ( helper_106 * helper_107 + helper_111 );
	result_0( 1, 0 ) = helper_87 * ( helper_104 + helper_116 * helper_99 );
	result_0( 1, 1 ) = helper_56 * ( -pow2( helper_117 ) * helper_83 + helper_117 * helper_65 * helper_92 + helper_117 * helper_84 * helper_91 + 3.0 );
	result_0( 1, 2 ) = helper_87 * ( -helper_105 * v1_x - helper_107 * helper_116 + helper_118 );
	result_0( 2, 0 ) = helper_87 * ( -helper_105 * v1_y + helper_111 + helper_119 * helper_99 );
	result_0( 2, 1 ) = helper_87 * ( helper_118 - helper_119 * helper_91 );
	result_0( 2, 2 ) = helper_56 * ( -helper_108 * helper_109 * helper_65 - 1.11111111111111 * pow2( helper_109 ) * helper_84 + 3.0 );
}

void floatTetWild::AMIPS_hessian( const std::array<Scalar, 12>& T, Matrix3& result_0 )
{
#if 0
	AMIPS_hessian_v2( T, result_0 );
#else
	Matrix3 matOld, matNew;
	AMIPS_hessian_v1( T, matOld );
	AMIPS_hessian_v2( T, matNew );
	Matrix3 diff = matNew - matOld;
	__debugbreak();
	result_0 = matOld;
#endif
}