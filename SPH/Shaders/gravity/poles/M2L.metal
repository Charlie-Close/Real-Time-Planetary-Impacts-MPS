//
//  M2L.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include <metal_stdlib>
#include "poles.h"
using namespace metal;

Local M2L(float3 x, Multipole mp) {
    float3 r = x - mp.pos;
    Derivatives deriv = derivatives(r, GRAVITY_SMOOTHING_LENGTH);
    Local local;
    
    const float M_000 = mp.expansion[M];
    const float D_000 = deriv.expansion[M];
    local.expansion[M] = M_000 * D_000;
    
#if P > 0
    const float D_100 = deriv.expansion[X];
    const float D_010 = deriv.expansion[Y];
    const float D_001 = deriv.expansion[Z];

    /*  1st order multipole term (addition to rank 1)*/
    local.expansion[X] = M_000 * D_100;
    local.expansion[Y] = M_000 * D_010;
    local.expansion[Z] = M_000 * D_001;
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
    local.expansion[XX] = M_000 * D_200;
    local.expansion[YY] = M_000 * D_020;
    local.expansion[ZZ] = M_000 * D_002;
    local.expansion[XY] = M_000 * D_110;
    local.expansion[XZ] = M_000 * D_101;
    local.expansion[YZ] = M_000 * D_011;
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
    local.expansion[XXX] = M_000 * D_300;
    local.expansion[YYY] = M_000 * D_030;
    local.expansion[ZZZ] = M_000 * D_003;
    local.expansion[XXY] = M_000 * D_210;
    local.expansion[XXZ] = M_000 * D_201;
    local.expansion[XYY] = M_000 * D_120;
    local.expansion[YYZ] = M_000 * D_021;
    local.expansion[XZZ] = M_000 * D_102;
    local.expansion[YZZ] = M_000 * D_012;
    local.expansion[XYZ] = M_000 * D_111;
#endif
    return local;
}


