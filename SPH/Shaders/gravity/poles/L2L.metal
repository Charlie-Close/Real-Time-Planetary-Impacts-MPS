//
//  L2L.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include "poles.h"
#include <metal_stdlib>
using namespace metal;

Local transformLocal(Local local, float3 r) {
    Local newLocal;
    newLocal.pos = local.pos + r;
    
    newLocal.expansion[M] = local.expansion[M];
    
    const float X_100 = r.x;
    const float X_010 = r.y;
    const float X_001 = r.z;
    
    newLocal.expansion[M] += X_100 * local.expansion[X] + X_010 * local.expansion[Y] + X_001 * local.expansion[Z];
    newLocal.expansion[X] = local.expansion[X];
    newLocal.expansion[Y] = local.expansion[Y];
    newLocal.expansion[Z] = local.expansion[Z];
    
#if P > 1
    const float R_200 = r.x * r.x;
    const float R_020 = r.y * r.y;
    const float R_002 = r.z * r.z;
    const float R_110 = r.x * r.y;
    const float R_101 = r.x * r.z;
    const float R_011 = r.y * r.z;
    
    const float X_200 = 0.5 * R_200;
    const float X_020 = 0.5 * R_020;
    const float X_002 = 0.5 * R_002;
    const float X_110 = R_110;
    const float X_101 = R_101;
    const float X_011 = R_011;
    
    newLocal.expansion[M] += X_200 * local.expansion[XX] + X_020 * local.expansion[YY] + X_002 * local.expansion[ZZ];
    newLocal.expansion[M] += X_110 * local.expansion[XY] + X_101 * local.expansion[XZ] + X_011 * local.expansion[YZ];
    newLocal.expansion[X] += X_100 * local.expansion[XX] + X_010 * local.expansion[XY] + X_001 * local.expansion[XZ];
    newLocal.expansion[Y] += X_100 * local.expansion[XY] + X_010 * local.expansion[YY] + X_001 * local.expansion[YZ];
    newLocal.expansion[Z] += X_100 * local.expansion[XZ] + X_010 * local.expansion[YZ] + X_001 * local.expansion[ZZ];
    
    newLocal.expansion[XX] = local.expansion[XX];
    newLocal.expansion[XY] = local.expansion[XY];
    newLocal.expansion[XZ] = local.expansion[XZ];
    newLocal.expansion[YY] = local.expansion[YY];
    newLocal.expansion[YZ] = local.expansion[YZ];
    newLocal.expansion[ZZ] = local.expansion[ZZ];
#endif
#if P > 2
    const float R_300 = r.x * R_200;
    const float R_030 = r.y * R_020;
    const float R_003 = r.z * R_002;
    const float R_210 = r.x * R_110;
    const float R_201 = r.x * R_101;
    const float R_120 = r.y * R_110;
    const float R_021 = r.y * R_011;
    const float R_102 = r.z * R_101;
    const float R_012 = r.z * R_011;
    const float R_111 = r.x * X_011;

    const float X_300 = 0.166666667 * R_300;
    const float X_030 = 0.166666667 * R_030;
    const float X_003 = 0.166666667 * R_003;
    const float X_210 = 0.5 * R_210;
    const float X_201 = 0.5 * R_201;
    const float X_120 = 0.5 * R_120;
    const float X_021 = 0.5 * R_021;
    const float X_102 = 0.5 * R_102;
    const float X_012 = 0.5 * R_012;
    const float X_111 = R_111;
    
    /* Shift 3rd order multipole term (addition to rank 0)*/
    newLocal.expansion[M] +=
        X_300 * local.expansion[XXX] + X_030 * local.expansion[YYY] + X_003 * local.expansion[ZZZ];
    newLocal.expansion[M] +=
        X_210 * local.expansion[XXY] + X_201 * local.expansion[XXZ] + X_120 * local.expansion[XYY];
    newLocal.expansion[M] +=
        X_021 * local.expansion[YYZ] + X_102 * local.expansion[XZZ] + X_012 * local.expansion[YZZ];
    newLocal.expansion[M] += X_111 * local.expansion[XYZ];

    /* Shift 3rd order multipole term (addition to rank 1)*/
    newLocal.expansion[X] +=
        X_200 * local.expansion[XXX] + X_020 * local.expansion[XYY] + X_002 * local.expansion[XZZ];
    newLocal.expansion[X] +=
        X_110 * local.expansion[XXY] + X_101 * local.expansion[XXZ] + X_011 * local.expansion[XYZ];
    newLocal.expansion[Y] +=
        X_200 * local.expansion[XXY] + X_020 * local.expansion[YYY] + X_002 * local.expansion[YZZ];
    newLocal.expansion[Y] +=
        X_110 * local.expansion[XYY] + X_101 * local.expansion[XYZ] + X_011 * local.expansion[YYZ];
    newLocal.expansion[Z] +=
        X_200 * local.expansion[XXZ] + X_020 * local.expansion[YYZ] + X_002 * local.expansion[ZZZ];
    newLocal.expansion[Z] +=
        X_110 * local.expansion[XYZ] + X_101 * local.expansion[XZZ] + X_011 * local.expansion[YZZ];

    /* Shift 3rd order multipole term (addition to rank 2)*/
    newLocal.expansion[XX] +=
        X_100 * local.expansion[XXX] + X_010 * local.expansion[XXY] + X_001 * local.expansion[XXZ];
    newLocal.expansion[YY] +=
        X_100 * local.expansion[XYY] + X_010 * local.expansion[YYY] + X_001 * local.expansion[YYZ];
    newLocal.expansion[ZZ] +=
        X_100 * local.expansion[XZZ] + X_010 * local.expansion[YZZ] + X_001 * local.expansion[ZZZ];
    newLocal.expansion[XY] +=
        X_100 * local.expansion[XXY] + X_010 * local.expansion[XYY] + X_001 * local.expansion[XYZ];
    newLocal.expansion[XZ] +=
        X_100 * local.expansion[XXZ] + X_010 * local.expansion[XYZ] + X_001 * local.expansion[XZZ];
    newLocal.expansion[YZ] +=
        X_100 * local.expansion[XYZ] + X_010 * local.expansion[YYZ] + X_001 * local.expansion[YZZ];

    /* Shift 3rd order multipole term (addition to rank 3)*/
    newLocal.expansion[XXX] = local.expansion[XXX];//F_300;
    newLocal.expansion[YYY] = local.expansion[YYY];//F_030;
    newLocal.expansion[ZZZ] = local.expansion[ZZZ];//F_003;
    newLocal.expansion[XXY] = local.expansion[XXY];//F_210;
    newLocal.expansion[XXZ] = local.expansion[XXZ];//F_201;
    newLocal.expansion[XYY] = local.expansion[XYY];//F_120;
    newLocal.expansion[YYZ] = local.expansion[YYZ];//F_021;
    newLocal.expansion[XZZ] = local.expansion[XZZ];//F_102;
    newLocal.expansion[YZZ] = local.expansion[YZZ];//F_012;
    newLocal.expansion[XYZ] = local.expansion[XYZ];//F_111;
#endif
#if P > 3
    const float R_400 = r.x * R_300;
    const float R_040 = r.y * R_030;
    const float R_004 = r.z * R_003;
    const float R_310 = r.x * R_210;
    const float R_301 = r.x * R_201;
    const float R_130 = r.y * R_120;
    const float R_031 = r.y * R_021;
    const float R_103 = r.z * R_102;
    const float R_013 = r.z * R_012;
    const float R_220 = R_110 * R_110;
    const float R_202 = R_101 * R_101;
    const float R_022 = R_011 * R_011;
    const float R_211 = r.x * R_111;
    const float R_121 = r.y * R_111;
    const float R_112 = r.z * R_111;
    
    const float X_400 = 0.041666666666666667 * R_400;
    const float X_040 = 0.041666666666666667 * R_040;
    const float X_004 = 0.041666666666666667 * R_004;
    const float X_310 = 0.1666666666666667 * R_310;
    const float X_301 = 0.1666666666666667 * R_301;
    const float X_130 = 0.1666666666666667 * R_130;
    const float X_031 = 0.1666666666666667 * R_031;
    const float X_103 = 0.1666666666666667 * R_103;
    const float X_013 = 0.1666666666666667 * R_013;
    const float X_220 = 0.25 * R_220;
    const float X_202 = 0.25 * R_202;
    const float X_022 = 0.25 * R_022;
    const float X_211 = 0.5 * R_211;
    const float X_121 = 0.5 * R_121;
    const float X_112 = 0.5 * R_112;

    newLocal.expansion[M] +=
        X_004 * local.expansion[ZZZZ] + X_013 * local.expansion[YZZZ] + X_022 * local.expansion[YYZZ] +
        X_031 * local.expansion[YYYZ] + X_040 * local.expansion[YYYY] + X_103 * local.expansion[XZZZ] +
        X_112 * local.expansion[XYZZ] + X_121 * local.expansion[XYYZ] + X_130 * local.expansion[XYYY] +
        X_202 * local.expansion[XXZZ] + X_211 * local.expansion[XXYZ] + X_220 * local.expansion[XXYY] +
        X_301 * local.expansion[XXXZ] + X_310 * local.expansion[XXXY] + X_400 * local.expansion[XXXX];

    /* Shift 4th order field tensor terms (addition to rank 1) */
    newLocal.expansion[Z] += X_003 * local.expansion[ZZZZ] + X_012 * local.expansion[YZZZ] +
                 X_021 * local.expansion[YYZZ] + X_030 * local.expansion[YYYZ] +
                 X_102 * local.expansion[XZZZ] + X_111 * local.expansion[XYZZ] +
                 X_120 * local.expansion[XYYZ] + X_201 * local.expansion[XXZZ] +
                 X_210 * local.expansion[XXYZ] + X_300 * local.expansion[XXXZ];
    newLocal.expansion[Y] += X_003 * local.expansion[YZZZ] + X_012 * local.expansion[YYZZ] +
                 X_021 * local.expansion[YYYZ] + X_030 * local.expansion[YYYY] +
                 X_102 * local.expansion[XYZZ] + X_111 * local.expansion[XYYZ] +
                 X_120 * local.expansion[XYYY] + X_201 * local.expansion[XXYZ] +
                 X_210 * local.expansion[XXYY] + X_300 * local.expansion[XXXY];
    newLocal.expansion[X] += X_003 * local.expansion[XZZZ] + X_012 * local.expansion[XYZZ] +
                 X_021 * local.expansion[XYYZ] + X_030 * local.expansion[XYYY] +
                 X_102 * local.expansion[XXZZ] + X_111 * local.expansion[XXYZ] +
                 X_120 * local.expansion[XXYY] + X_201 * local.expansion[XXXZ] +
                 X_210 * local.expansion[XXXY] + X_300 * local.expansion[XXXX];

    /* Shift 4th order field tensor terms (addition to rank 2) */
    newLocal.expansion[ZZ] += X_002 * local.expansion[ZZZZ] + X_011 * local.expansion[YZZZ] +
                 X_020 * local.expansion[YYZZ] + X_101 * local.expansion[XZZZ] +
                 X_110 * local.expansion[XYZZ] + X_200 * local.expansion[XXZZ];
    newLocal.expansion[YZ] += X_002 * local.expansion[YZZZ] + X_011 * local.expansion[YYZZ] +
                 X_020 * local.expansion[YYYZ] + X_101 * local.expansion[XYZZ] +
                 X_110 * local.expansion[XYYZ] + X_200 * local.expansion[XXYZ];
    newLocal.expansion[YY] += X_002 * local.expansion[YYZZ] + X_011 * local.expansion[YYYZ] +
                 X_020 * local.expansion[YYYY] + X_101 * local.expansion[XYYZ] +
                 X_110 * local.expansion[XYYY] + X_200 * local.expansion[XXYY];
    newLocal.expansion[XZ] += X_002 * local.expansion[XZZZ] + X_011 * local.expansion[XYZZ] +
                 X_020 * local.expansion[XYYZ] + X_101 * local.expansion[XXZZ] +
                 X_110 * local.expansion[XXYZ] + X_200 * local.expansion[XXXZ];
    newLocal.expansion[XY] += X_002 * local.expansion[XYZZ] + X_011 * local.expansion[XYYZ] +
                 X_020 * local.expansion[XYYY] + X_101 * local.expansion[XXYZ] +
                 X_110 * local.expansion[XXYY] + X_200 * local.expansion[XXXY];
    newLocal.expansion[XX] += X_002 * local.expansion[XXZZ] + X_011 * local.expansion[XXYZ] +
                 X_020 * local.expansion[XXYY] + X_101 * local.expansion[XXXZ] +
                 X_110 * local.expansion[XXXY] + X_200 * local.expansion[XXXX];

    /* Shift 4th order field tensor terms (addition to rank 3) */
    newLocal.expansion[ZZZ] +=
        X_001 * local.expansion[ZZZZ] + X_010 * local.expansion[YZZZ] + X_100 * local.expansion[XZZZ];
    newLocal.expansion[YZZ] +=
        X_001 * local.expansion[YZZZ] + X_010 * local.expansion[YYZZ] + X_100 * local.expansion[XYZZ];
    newLocal.expansion[YYZ] +=
        X_001 * local.expansion[YYZZ] + X_010 * local.expansion[YYYZ] + X_100 * local.expansion[XYYZ];
    newLocal.expansion[YYY] +=
        X_001 * local.expansion[YYYZ] + X_010 * local.expansion[YYYY] + X_100 * local.expansion[XYYY];
    newLocal.expansion[XZZ] +=
        X_001 * local.expansion[XZZZ] + X_010 * local.expansion[XYZZ] + X_100 * local.expansion[XXZZ];
    newLocal.expansion[XYZ] +=
        X_001 * local.expansion[XYZZ] + X_010 * local.expansion[XYYZ] + X_100 * local.expansion[XXYZ];
    newLocal.expansion[XYY] +=
        X_001 * local.expansion[XYYZ] + X_010 * local.expansion[XYYY] + X_100 * local.expansion[XXYY];
    newLocal.expansion[XXZ] +=
        X_001 * local.expansion[XXZZ] + X_010 * local.expansion[XXYZ] + X_100 * local.expansion[XXXZ];
    newLocal.expansion[XXY] +=
        X_001 * local.expansion[XXYZ] + X_010 * local.expansion[XXYY] + X_100 * local.expansion[XXXY];
    newLocal.expansion[XXX] +=
        X_001 * local.expansion[XXXZ] + X_010 * local.expansion[XXXY] + X_100 * local.expansion[XXXX];

    /* Shift 4th order field tensor terms (addition to rank 4) */
    newLocal.expansion[ZZZZ] = local.expansion[ZZZZ];
    newLocal.expansion[YZZZ] = local.expansion[YZZZ];
    newLocal.expansion[YYZZ] = local.expansion[YYZZ];
    newLocal.expansion[YYYZ] = local.expansion[YYYZ];
    newLocal.expansion[YYYY] = local.expansion[YYYY];
    newLocal.expansion[XZZZ] = local.expansion[XZZZ];
    newLocal.expansion[XYZZ] = local.expansion[XYZZ];
    newLocal.expansion[XYYZ] = local.expansion[XYYZ];
    newLocal.expansion[XYYY] = local.expansion[XYYY];
    newLocal.expansion[XXZZ] = local.expansion[XXZZ];
    newLocal.expansion[XXYZ] = local.expansion[XXYZ];
    newLocal.expansion[XXYY] = local.expansion[XXYY];
    newLocal.expansion[XXXZ] = local.expansion[XXXZ];
    newLocal.expansion[XXXY] = local.expansion[XXXY];
    newLocal.expansion[XXXX] = local.expansion[XXXX];
#endif
    return newLocal;
}

void L2L(device int* treeStructure, device Multipole* multipoles, device Local* locals, thread Multipole& mp, thread Local& local, int treePointer) {
    // Give our local expansion to our children
    int start = treePointer + 2;
    int end = start + 8;
    for (int j = start; j < end; j++) {
        int nodePointer = treeStructure[j];
        int nodeDataPointer = treeStructure[nodePointer + 1];
        Multipole nodeMp = multipoles[nodeDataPointer];
        float3 r = nodeMp.pos - local.pos;
        locals[nodeDataPointer] = transformLocal(local, r);
    }
}
