//
//  L2P.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include "poles.h"
#include <metal_stdlib>
using namespace metal;

float3 L2P(Local local, float3 pos) {
    float3 acc = { 0, 0, 0 };
    float3 r = pos - local.pos;
    
    acc += float3(local.expansion[X], local.expansion[Y], local.expansion[Z]);
#if P > 1
    acc.x += r.x * local.expansion[XX] + r.y * local.expansion[XY] + r.z * local.expansion[XZ];
    acc.y += r.x * local.expansion[XY] + r.y * local.expansion[YY] + r.z * local.expansion[YZ];
    acc.z += r.x * local.expansion[XZ] + r.y * local.expansion[YZ] + r.z * local.expansion[ZZ];
#endif
#if P > 2
    acc.x += 0.5 * r.x * r.x * local.expansion[XXX] + 0.5 * r.y * r.y * local.expansion[XYY] + 0.5 * r.z * r.z * local.expansion[XZZ];
    acc.x += r.x * r.y * local.expansion[XXY] + r.x * r.z * local.expansion[XXZ] + r.y * r.z * local.expansion[XYZ];
    acc.y += 0.5 * r.x * r.x * local.expansion[XXY] + 0.5 * r.y * r.y * local.expansion[YYY] + 0.5 * r.z * r.z * local.expansion[YZZ];
    acc.y += r.x * r.y * local.expansion[XYY] + r.x * r.z * local.expansion[XYZ] + r.y * r.z * local.expansion[YYZ];
    acc.z += 0.5 * r.x * r.x * local.expansion[XXZ] + 0.5 * r.y * r.y * local.expansion[YYZ] + 0.5 * r.z * r.z * local.expansion[ZZZ];
    acc.z += r.x * r.y * local.expansion[XYZ] + r.x * r.z * local.expansion[XZZ] + r.y * r.z * local.expansion[YZZ];
#endif
 
    return acc;
}
