//
//  M2L.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include <metal_stdlib>
#include "poles.h"
using namespace metal;

void M2L(thread Multipole& mp, thread Local& local) {
    float3 r = local.pos - mp.pos;
    Derivatives deriv = derivatives(r, mp.eta);
    
    const float M_000 = mp.expansion[M];
    const float D_000 = deriv.expansion[M];
    
    local.expansion[M] += M_000 * D_000;
    
#if P > 0
    const float D_100 = deriv.expansion[X];
    const float D_010 = deriv.expansion[Y];
    const float D_001 = deriv.expansion[Z];

    /*  1st order multipole term (addition to rank 1)*/
    local.expansion[X] += M_000 * D_100;
    local.expansion[Y] += M_000 * D_010;
    local.expansion[Z] += M_000 * D_001;
#endif
#if P > 1
    const float M_200 = mp.expansion[XX];
    const float M_020 = mp.expansion[YY];
    const float M_002 = mp.expansion[ZZ];
    const float M_110 = mp.expansion[XY];
    const float M_101 = mp.expansion[XZ];
    const float M_011 = mp.expansion[YZ];//M_011;

    const float D_200 = deriv.expansion[XX];
    const float D_020 = deriv.expansion[YY];
    const float D_002 = deriv.expansion[ZZ];
    const float D_110 = deriv.expansion[XY];
    const float D_101 = deriv.expansion[XZ];
    const float D_011 = deriv.expansion[YZ];

    /*  2nd order multipole term (addition to rank 0)*/
    local.expansion[M] += M_200 * D_200 + M_020 * D_020 + M_002 * D_002;
    local.expansion[M] += M_110 * D_110 + M_101 * D_101 + M_011 * D_011;

    /*  2nd order multipole term (addition to rank 2)*/
    local.expansion[XX] += M_000 * D_200;
    local.expansion[YY] += M_000 * D_020;
    local.expansion[ZZ] += M_000 * D_002;
    local.expansion[XY] += M_000 * D_110;
    local.expansion[XZ] += M_000 * D_101;
    local.expansion[YZ] += M_000 * D_011;
#endif
#if P > 2
    const float M_300 = mp.expansion[XXX];//M_300;
    const float M_030 = mp.expansion[YYY];//M_030;
    const float M_003 = mp.expansion[ZZZ];//M_003;
    const float M_210 = mp.expansion[XXY];//M_210;
    const float M_201 = mp.expansion[XXZ];//M_201;
    const float M_021 = mp.expansion[YYZ];//M_021;
    const float M_120 = mp.expansion[XYY];//M_120;
    const float M_012 = mp.expansion[YZZ];//M_012;
    const float M_102 = mp.expansion[XZZ];//M_102;
    const float M_111 = mp.expansion[XYZ];//M_111;

    const float D_300 = deriv.expansion[XXX];//D_300;
    const float D_030 = deriv.expansion[YYY];//D_030;
    const float D_003 = deriv.expansion[ZZZ];//D_003;
    const float D_210 = deriv.expansion[XXY];//D_210;
    const float D_201 = deriv.expansion[XXZ];//D_201;
    const float D_021 = deriv.expansion[YYZ];//D_021;
    const float D_120 = deriv.expansion[XYY];//D_120;
    const float D_012 = deriv.expansion[YZZ];//D_012;
    const float D_102 = deriv.expansion[XZZ];//D_102;
    const float D_111 = deriv.expansion[XYZ];//D_111;

    /*  3rd order multipole term (addition to rank 0)*/
    local.expansion[M] += M_300 * D_300 + M_030 * D_030 + M_003 * D_003;
    local.expansion[M] += M_210 * D_210 + M_201 * D_201 + M_120 * D_120;
    local.expansion[M] += M_021 * D_021 + M_102 * D_102 + M_012 * D_012;
    local.expansion[M] += M_111 * D_111;

    /*  3rd order multipole term (addition to rank 1)*/
    local.expansion[X] += M_200 * D_300 + M_020 * D_120 + M_002 * D_102;
    local.expansion[X] += M_110 * D_210 + M_101 * D_201 + M_011 * D_111;
    local.expansion[Y] += M_200 * D_210 + M_020 * D_030 + M_002 * D_012;
    local.expansion[Y] += M_110 * D_120 + M_101 * D_111 + M_011 * D_021;
    local.expansion[Z] += M_200 * D_201 + M_020 * D_021 + M_002 * D_003;
    local.expansion[Z] += M_110 * D_111 + M_101 * D_102 + M_011 * D_012;

    /*  3rd order multipole term (addition to rank 3)*/
    local.expansion[XXX] += M_000 * D_300;
    local.expansion[YYY] += M_000 * D_030;
    local.expansion[ZZZ] += M_000 * D_003;
    local.expansion[XXY] += M_000 * D_210;
    local.expansion[XXZ] += M_000 * D_201;
    local.expansion[XYY] += M_000 * D_120;
    local.expansion[YYZ] += M_000 * D_021;
    local.expansion[XZZ] += M_000 * D_102;
    local.expansion[YZZ] += M_000 * D_012;
    local.expansion[XYZ] += M_000 * D_111;
#endif
#if P > 3
    const float M_400 = mp.expansion[XXXX];//M_400;
    const float M_040 = mp.expansion[YYYY];//M_040;
    const float M_004 = mp.expansion[ZZZZ];//M_004;
    const float M_310 = mp.expansion[XXXY];//M_310;
    const float M_301 = mp.expansion[XXXZ];//M_301;
    const float M_031 = mp.expansion[YYYZ];//M_031;
    const float M_130 = mp.expansion[XYYY];//M_130;
    const float M_013 = mp.expansion[YZZZ];//M_013;
    const float M_103 = mp.expansion[XZZZ];//M_103;
    const float M_220 = mp.expansion[XXYY];//M_220;
    const float M_202 = mp.expansion[XXZZ];//M_202;
    const float M_022 = mp.expansion[YYZZ];//M_022;
    const float M_211 = mp.expansion[XXYZ];//M_211;
    const float M_121 = mp.expansion[XYYZ];//M_121;
    const float M_112 = mp.expansion[XYZZ];//M_112;

    const float D_400 = deriv.expansion[XXXX];//D_400;
    const float D_040 = deriv.expansion[YYYY];//D_040;
    const float D_004 = deriv.expansion[ZZZZ];//D_004;
    const float D_310 = deriv.expansion[XXXY];//D_310;
    const float D_301 = deriv.expansion[XXXZ];//D_301;
    const float D_031 = deriv.expansion[YYYZ];//D_031;
    const float D_130 = deriv.expansion[XYYY];//D_130;
    const float D_013 = deriv.expansion[YZZZ];//D_013;
    const float D_103 = deriv.expansion[XZZZ];//D_103;
    const float D_220 = deriv.expansion[XXYY];//D_220;
    const float D_202 = deriv.expansion[XXZZ];//D_202;
    const float D_022 = deriv.expansion[YYZZ];//D_022;
    const float D_211 = deriv.expansion[XXYZ];//D_211;
    const float D_121 = deriv.expansion[XYYZ];//D_121;
    const float D_112 = deriv.expansion[XYZZ];//D_112;
    /* Compute 4th order field tensor terms (addition to rank 0) */
    local.expansion[M] += M_004 * D_004 + M_013 * D_013 + M_022 * D_022 + M_031 * D_031 +
                M_040 * D_040 + M_103 * D_103 + M_112 * D_112 + M_121 * D_121 +
                M_130 * D_130 + M_202 * D_202 + M_211 * D_211 + M_220 * D_220 +
                M_301 * D_301 + M_310 * D_310 + M_400 * D_400;

    /* Compute 4th order field tensor terms (addition to rank 1) */
    local.expansion[Z] += M_003 * D_004 + M_012 * D_013 + M_021 * D_022 + M_030 * D_031 +
                M_102 * D_103 + M_111 * D_112 + M_120 * D_121 + M_201 * D_202 +
                M_210 * D_211 + M_300 * D_301;
    local.expansion[Y] += M_003 * D_013 + M_012 * D_022 + M_021 * D_031 + M_030 * D_040 +
                M_102 * D_112 + M_111 * D_121 + M_120 * D_130 + M_201 * D_211 +
                M_210 * D_220 + M_300 * D_310;
    local.expansion[X] += M_003 * D_103 + M_012 * D_112 + M_021 * D_121 + M_030 * D_130 +
                M_102 * D_202 + M_111 * D_211 + M_120 * D_220 + M_201 * D_301 +
                M_210 * D_310 + M_300 * D_400;

    /* Compute 4th order field tensor terms (addition to rank 2) */
    local.expansion[ZZ] += M_002 * D_004 + M_011 * D_013 + M_020 * D_022 + M_101 * D_103 +
                M_110 * D_112 + M_200 * D_202;
    local.expansion[YZ] += M_002 * D_013 + M_011 * D_022 + M_020 * D_031 + M_101 * D_112 +
                M_110 * D_121 + M_200 * D_211;
    local.expansion[YY] += M_002 * D_022 + M_011 * D_031 + M_020 * D_040 + M_101 * D_121 +
                M_110 * D_130 + M_200 * D_220;
    local.expansion[XZ] += M_002 * D_103 + M_011 * D_112 + M_020 * D_121 + M_101 * D_202 +
                M_110 * D_211 + M_200 * D_301;
    local.expansion[XY] += M_002 * D_112 + M_011 * D_121 + M_020 * D_130 + M_101 * D_211 +
                M_110 * D_220 + M_200 * D_310;
    local.expansion[XX] += M_002 * D_202 + M_011 * D_211 + M_020 * D_220 + M_101 * D_301 +
                M_110 * D_310 + M_200 * D_400;

    /* Compute 4th order field tensor terms (addition to rank 4) */
    local.expansion[ZZZZ] += M_000 * D_004;
    local.expansion[YZZZ] += M_000 * D_013;
    local.expansion[YYZZ] += M_000 * D_022;
    local.expansion[YYYZ] += M_000 * D_031;
    local.expansion[YYYY] += M_000 * D_040;
    local.expansion[XZZZ] += M_000 * D_103;
    local.expansion[XYZZ] += M_000 * D_112;
    local.expansion[XYYZ] += M_000 * D_121;
    local.expansion[XYYY] += M_000 * D_130;
    local.expansion[XXZZ] += M_000 * D_202;
    local.expansion[XXYZ] += M_000 * D_211;
    local.expansion[XXYY] += M_000 * D_220;
    local.expansion[XXXZ] += M_000 * D_301;
    local.expansion[XXXY] += M_000 * D_310;
    local.expansion[XXXX] += M_000 * D_400;
#endif
}


//Local M2L(float3 x, Multipole mp) {
//    float3 r = x - mp.pos;
//    Derivatives deriv = derivatives(r, mp.eta);
//    Local local;
//    
//    const float M_000 = mp.expansion[M];
//    
//    const float D_000 = deriv.expansion[M];
//    local.expansion[M] = M_000 * D_000;
//    
//#if P > 0
//    const float D_100 = deriv.expansion[X];
//    const float D_010 = deriv.expansion[Y];
//    const float D_001 = deriv.expansion[Z];
//
//    /*  1st order multipole term (addition to rank 1)*/
//    local.expansion[X] = M_000 * D_100;
//    local.expansion[Y] = M_000 * D_010;
//    local.expansion[Z] = M_000 * D_001;
//#endif
//#if P > 1
//    const float M_200 = mp.expansion[XX];
//    const float M_020 = mp.expansion[YY];
//    const float M_002 = mp.expansion[ZZ];
//    const float M_110 = mp.expansion[XY];
//    const float M_101 = mp.expansion[XZ];
//    const float M_011 = mp.expansion[YZ];//M_011;
//
//    const float D_200 = deriv.expansion[XX];
//    const float D_020 = deriv.expansion[YY];
//    const float D_002 = deriv.expansion[ZZ];
//    const float D_110 = deriv.expansion[XY];
//    const float D_101 = deriv.expansion[XZ];
//    const float D_011 = deriv.expansion[YZ];
//
//    /*  2nd order multipole term (addition to rank 0)*/
//    local.expansion[M] += M_200 * D_200 + M_020 * D_020 + M_002 * D_002;
//    local.expansion[M] += M_110 * D_110 + M_101 * D_101 + M_011 * D_011;
//
//    /*  2nd order multipole term (addition to rank 2)*/
//    local.expansion[XX] = M_000 * D_200;
//    local.expansion[YY] = M_000 * D_020;
//    local.expansion[ZZ] = M_000 * D_002;
//    local.expansion[XY] = M_000 * D_110;
//    local.expansion[XZ] = M_000 * D_101;
//    local.expansion[YZ] = M_000 * D_011;
//#endif
//#if P > 2
//    const float M_300 = mp.expansion[XXX];//M_300;
//    const float M_030 = mp.expansion[YYY];//M_030;
//    const float M_003 = mp.expansion[ZZZ];//M_003;
//    const float M_210 = mp.expansion[XXY];//M_210;
//    const float M_201 = mp.expansion[XXZ];//M_201;
//    const float M_021 = mp.expansion[YYZ];//M_021;
//    const float M_120 = mp.expansion[XYY];//M_120;
//    const float M_012 = mp.expansion[YZZ];//M_012;
//    const float M_102 = mp.expansion[XZZ];//M_102;
//    const float M_111 = mp.expansion[XYZ];//M_111;
//
//    const float D_300 = deriv.expansion[XXX];//D_300;
//    const float D_030 = deriv.expansion[YYY];//D_030;
//    const float D_003 = deriv.expansion[ZZZ];//D_003;
//    const float D_210 = deriv.expansion[XXY];//D_210;
//    const float D_201 = deriv.expansion[XXZ];//D_201;
//    const float D_021 = deriv.expansion[YYZ];//D_021;
//    const float D_120 = deriv.expansion[XYY];//D_120;
//    const float D_012 = deriv.expansion[YZZ];//D_012;
//    const float D_102 = deriv.expansion[XZZ];//D_102;
//    const float D_111 = deriv.expansion[XYZ];//D_111;
//
//    /*  3rd order multipole term (addition to rank 0)*/
//    local.expansion[M] += M_300 * D_300 + M_030 * D_030 + M_003 * D_003;
//    local.expansion[M] += M_210 * D_210 + M_201 * D_201 + M_120 * D_120;
//    local.expansion[M] += M_021 * D_021 + M_102 * D_102 + M_012 * D_012;
//    local.expansion[M] += M_111 * D_111;
//
//    /*  3rd order multipole term (addition to rank 1)*/
//    local.expansion[X] += M_200 * D_300 + M_020 * D_120 + M_002 * D_102;
//    local.expansion[X] += M_110 * D_210 + M_101 * D_201 + M_011 * D_111;
//    local.expansion[Y] += M_200 * D_210 + M_020 * D_030 + M_002 * D_012;
//    local.expansion[Y] += M_110 * D_120 + M_101 * D_111 + M_011 * D_021;
//    local.expansion[Z] += M_200 * D_201 + M_020 * D_021 + M_002 * D_003;
//    local.expansion[Z] += M_110 * D_111 + M_101 * D_102 + M_011 * D_012;
//
//    /*  3rd order multipole term (addition to rank 3)*/
//    local.expansion[XXX] = M_000 * D_300;
//    local.expansion[YYY] = M_000 * D_030;
//    local.expansion[ZZZ] = M_000 * D_003;
//    local.expansion[XXY] = M_000 * D_210;
//    local.expansion[XXZ] = M_000 * D_201;
//    local.expansion[XYY] = M_000 * D_120;
//    local.expansion[YYZ] = M_000 * D_021;
//    local.expansion[XZZ] = M_000 * D_102;
//    local.expansion[YZZ] = M_000 * D_012;
//    local.expansion[XYZ] = M_000 * D_111;
//#endif
//#if P > 3
//    const float M_400 = mp.expansion[XXXX];//M_400;
//    const float M_040 = mp.expansion[YYYY];//M_040;
//    const float M_004 = mp.expansion[ZZZZ];//M_004;
//    const float M_310 = mp.expansion[XXXY];//M_310;
//    const float M_301 = mp.expansion[XXXZ];//M_301;
//    const float M_031 = mp.expansion[YYYZ];//M_031;
//    const float M_130 = mp.expansion[XYYY];//M_130;
//    const float M_013 = mp.expansion[YZZZ];//M_013;
//    const float M_103 = mp.expansion[XZZZ];//M_103;
//    const float M_220 = mp.expansion[XXYY];//M_220;
//    const float M_202 = mp.expansion[XXZZ];//M_202;
//    const float M_022 = mp.expansion[YYZZ];//M_022;
//    const float M_211 = mp.expansion[XXYZ];//M_211;
//    const float M_121 = mp.expansion[XYYZ];//M_121;
//    const float M_112 = mp.expansion[XYZZ];//M_112;
//
//    const float D_400 = deriv.expansion[XXXX];//D_400;
//    const float D_040 = deriv.expansion[YYYY];//D_040;
//    const float D_004 = deriv.expansion[ZZZZ];//D_004;
//    const float D_310 = deriv.expansion[XXXY];//D_310;
//    const float D_301 = deriv.expansion[XXXZ];//D_301;
//    const float D_031 = deriv.expansion[YYYZ];//D_031;
//    const float D_130 = deriv.expansion[XYYY];//D_130;
//    const float D_013 = deriv.expansion[YZZZ];//D_013;
//    const float D_103 = deriv.expansion[XZZZ];//D_103;
//    const float D_220 = deriv.expansion[XXYY];//D_220;
//    const float D_202 = deriv.expansion[XXZZ];//D_202;
//    const float D_022 = deriv.expansion[YYZZ];//D_022;
//    const float D_211 = deriv.expansion[XXYZ];//D_211;
//    const float D_121 = deriv.expansion[XYYZ];//D_121;
//    const float D_112 = deriv.expansion[XYZZ];//D_112;
//    /* Compute 4th order field tensor terms (addition to rank 0) */
//    local.expansion[M] += M_004 * D_004 + M_013 * D_013 + M_022 * D_022 + M_031 * D_031 +
//                M_040 * D_040 + M_103 * D_103 + M_112 * D_112 + M_121 * D_121 +
//                M_130 * D_130 + M_202 * D_202 + M_211 * D_211 + M_220 * D_220 +
//                M_301 * D_301 + M_310 * D_310 + M_400 * D_400;
//
//    /* Compute 4th order field tensor terms (addition to rank 1) */
//    local.expansion[Z] += M_003 * D_004 + M_012 * D_013 + M_021 * D_022 + M_030 * D_031 +
//                M_102 * D_103 + M_111 * D_112 + M_120 * D_121 + M_201 * D_202 +
//                M_210 * D_211 + M_300 * D_301;
//    local.expansion[Y] += M_003 * D_013 + M_012 * D_022 + M_021 * D_031 + M_030 * D_040 +
//                M_102 * D_112 + M_111 * D_121 + M_120 * D_130 + M_201 * D_211 +
//                M_210 * D_220 + M_300 * D_310;
//    local.expansion[X] += M_003 * D_103 + M_012 * D_112 + M_021 * D_121 + M_030 * D_130 +
//                M_102 * D_202 + M_111 * D_211 + M_120 * D_220 + M_201 * D_301 +
//                M_210 * D_310 + M_300 * D_400;
//
//    /* Compute 4th order field tensor terms (addition to rank 2) */
//    local.expansion[ZZ] += M_002 * D_004 + M_011 * D_013 + M_020 * D_022 + M_101 * D_103 +
//                M_110 * D_112 + M_200 * D_202;
//    local.expansion[YZ] += M_002 * D_013 + M_011 * D_022 + M_020 * D_031 + M_101 * D_112 +
//                M_110 * D_121 + M_200 * D_211;
//    local.expansion[YY] += M_002 * D_022 + M_011 * D_031 + M_020 * D_040 + M_101 * D_121 +
//                M_110 * D_130 + M_200 * D_220;
//    local.expansion[XZ] += M_002 * D_103 + M_011 * D_112 + M_020 * D_121 + M_101 * D_202 +
//                M_110 * D_211 + M_200 * D_301;
//    local.expansion[XY] += M_002 * D_112 + M_011 * D_121 + M_020 * D_130 + M_101 * D_211 +
//                M_110 * D_220 + M_200 * D_310;
//    local.expansion[XX] += M_002 * D_202 + M_011 * D_211 + M_020 * D_220 + M_101 * D_301 +
//                M_110 * D_310 + M_200 * D_400;
//
//    /* Compute 4th order field tensor terms (addition to rank 4) */
//    local.expansion[ZZZZ] = M_000 * D_004;
//    local.expansion[YZZZ] = M_000 * D_013;
//    local.expansion[YYZZ] = M_000 * D_022;
//    local.expansion[YYYZ] = M_000 * D_031;
//    local.expansion[YYYY] = M_000 * D_040;
//    local.expansion[XZZZ] = M_000 * D_103;
//    local.expansion[XYZZ] = M_000 * D_112;
//    local.expansion[XYYZ] = M_000 * D_121;
//    local.expansion[XYYY] = M_000 * D_130;
//    local.expansion[XXZZ] = M_000 * D_202;
//    local.expansion[XXYZ] = M_000 * D_211;
//    local.expansion[XXYY] = M_000 * D_220;
//    local.expansion[XXXZ] = M_000 * D_301;
//    local.expansion[XXXY] = M_000 * D_310;
//    local.expansion[XXXX] = M_000 * D_400;
//#endif
//    return local;
//}


