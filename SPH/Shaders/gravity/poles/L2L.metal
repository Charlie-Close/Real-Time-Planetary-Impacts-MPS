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
    
    newLocal.expansion[M] += r.x * local.expansion[X] + r.y * local.expansion[Y] + r.z * local.expansion[Z];
    newLocal.expansion[X] = local.expansion[X];
    newLocal.expansion[Y] = local.expansion[Y];
    newLocal.expansion[Z] = local.expansion[Z];
    
#if P > 1
    newLocal.expansion[M] += 0.5 * (r.x * r.x * local.expansion[XX] + r.y * r.y * local.expansion[YY] + r.z * r.z * local.expansion[ZZ]);
    newLocal.expansion[M] += r.x * r.y * local.expansion[XY] + r.x * r.z * local.expansion[XZ] + r.y * r.z * local.expansion[YZ];
    newLocal.expansion[X] += r.x * local.expansion[XX] + r.y * local.expansion[XY] + r.z * local.expansion[XZ];
    newLocal.expansion[Y] += r.x * local.expansion[XY] + r.y * local.expansion[YY] + r.z * local.expansion[YZ];
    newLocal.expansion[Z] += r.x * local.expansion[XZ] + r.y * local.expansion[YZ] + r.z * local.expansion[ZZ];
    
    newLocal.expansion[XX] = local.expansion[XX];
    newLocal.expansion[XY] = local.expansion[XY];
    newLocal.expansion[XZ] = local.expansion[XZ];
    newLocal.expansion[YY] = local.expansion[YY];
    newLocal.expansion[YZ] = local.expansion[YZ];
    newLocal.expansion[ZZ] = local.expansion[ZZ];
#endif
#if P > 2
    /* Shift 3rd order multipole term (addition to rank 0)*/
    newLocal.expansion[M] +=
        0.166666667 * r.x * r.x * r.x * local.expansion[XXX] + 0.166666667 * r.y * r.y * r.y * local.expansion[YYY] + 0.166666667 * r.z * r.z * r.z * local.expansion[ZZZ];
    newLocal.expansion[M] +=
        0.5 * r.x * r.x * r.y * local.expansion[XXY] + 0.5 * r.x * r.x * r.z * local.expansion[XXZ] + 0.5 * r.x * r.y * r.y * local.expansion[XYY];
    newLocal.expansion[M] +=
        0.5 * r.y * r.y * r.z * local.expansion[YYZ] + 0.5 * r.x * r.z * r.z * local.expansion[XZZ] + 0.5 * r.y * r.z * r.z * local.expansion[YZZ];
    newLocal.expansion[M] += r.x * r.y * r.z * local.expansion[XYZ];

    /* Shift 3rd order multipole term (addition to rank 1)*/
    newLocal.expansion[X] +=
        0.5 * r.x * r.x * local.expansion[XXX] + 0.5 * r.y * r.y * local.expansion[XYY] + 0.5 * r.z * r.z * local.expansion[XZZ];
    newLocal.expansion[X] +=
        r.x * r.y * local.expansion[XXY] + r.x * r.z * local.expansion[XXZ] + r.y * r.z * local.expansion[XYZ];
    newLocal.expansion[Y] +=
        0.5 * r.x * r.x * local.expansion[XXY] + 0.5 * r.y * r.y * local.expansion[YYY] + 0.5 * r.z * r.z * local.expansion[YZZ];
    newLocal.expansion[Y] +=
        r.x * r.y * local.expansion[XYY] + r.x * r.z * local.expansion[XYZ] + r.y * r.z * local.expansion[YYZ];
    newLocal.expansion[Z] +=
        0.5 * r.x * r.x * local.expansion[XXZ] + 0.5 * r.y * r.y * local.expansion[YYZ] + 0.5 * r.z * r.z * local.expansion[ZZZ];
    newLocal.expansion[Z] +=
        r.x * r.y * local.expansion[XYZ] + r.x * r.z * local.expansion[XZZ] + r.y * r.z * local.expansion[YZZ];

    /* Shift 3rd order multipole term (addition to rank 2)*/
    newLocal.expansion[XX] +=
        r.x * local.expansion[XXX] + r.y * local.expansion[XXY] + r.z * local.expansion[XXZ];
    newLocal.expansion[YY] +=
        r.x * local.expansion[XYY] + r.y * local.expansion[YYY] + r.z * local.expansion[YYZ];
    newLocal.expansion[ZZ] +=
        r.x * local.expansion[XZZ] + r.y * local.expansion[YZZ] + r.z * local.expansion[ZZZ];
    newLocal.expansion[XY] +=
        r.x * local.expansion[XXY] + r.y * local.expansion[XYY] + r.z * local.expansion[XYZ];
    newLocal.expansion[XZ] +=
        r.x * local.expansion[XXZ] + r.y * local.expansion[XYZ] + r.z * local.expansion[XZZ];
    newLocal.expansion[YZ] +=
        r.x * local.expansion[XYZ] + r.y * local.expansion[YYZ] + r.z * local.expansion[YZZ];

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
        if (!nodeMp.active) {
            continue;
        }
        float3 r = nodeMp.pos - local.pos;
        locals[nodeDataPointer] = transformLocal(local, r);
    }
}
