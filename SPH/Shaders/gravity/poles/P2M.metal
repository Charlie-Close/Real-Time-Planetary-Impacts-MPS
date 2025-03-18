//
//  P2M.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include "poles.h"
#include <metal_stdlib>
using namespace metal;

Multipole P2M(device int* treeStructure, device float* masses, device float3* positions, device float* grav, device float* h, device bool* active, device int* nextActiveTime, device int* globalTime, device int& dt, int treePointer) {
    Multipole mp;
    int nParticles = treeStructure[treePointer];
    int start = treePointer + 2;
    int end = start + nParticles;
    
    mp.pos = { 0, 0, 0 };
    for (uint i = 0; i < N_EXPANSION_TERMS; i++) {
        mp.expansion[i] = 0;
    }
    mp.minGrav = MAXFLOAT;
    mp.active = false;
    mp.max = float3(-MAXFLOAT);
    mp.min = float3(MAXFLOAT);
    mp.eta = 0;
    
    for (int i = start; i < end; i++) {
        int p = treeStructure[i];
        float p_mass = masses[p];
        float3 p_position = positions[p];
        float p_h = h[p];
        float p_eta = min(p_h * GAMMA * PLUMBER_EQUIVALENT, GRAVITY_SMOOTHING_LENGTH);
        int isActive = active[p];
        if (!mp.active and isActive) {
            mp.active = true;
        }
        mp.max = max(mp.max, p_position);
        mp.min = min(mp.min, p_position);
        mp.eta = max(mp.eta, p_eta);
        
        mp.expansion[M] += p_mass;
        mp.pos += p_position * p_mass;
        mp.minGrav = min(mp.minGrav, grav[p]);
    }
    
    mp.pos /= mp.expansion[M];
    
    for (int i = start; i < end; i++) {
        int p = treeStructure[i];
        float p_mass = masses[p];
        float3 p_position = positions[p];
        float3 r = p_position - mp.pos;
        
        mp.expansion[X] -= r.x * p_mass;
        mp.expansion[Y] -= r.y * p_mass;
        mp.expansion[Z] -= r.z * p_mass;
        
#if P > 1
        mp.expansion[XX] += p_mass * r.x * r.x;
        mp.expansion[XY] += p_mass * r.x * r.y;
        mp.expansion[XZ] += p_mass * r.x * r.z;
        mp.expansion[YY] += p_mass * r.y * r.y;
        mp.expansion[YZ] += p_mass * r.y * r.z;
        mp.expansion[ZZ] += p_mass * r.z * r.z;
#endif
#if P > 2
        mp.expansion[XXX] -= p_mass * r.x * r.x * r.x;
        mp.expansion[XXY] -= p_mass * r.x * r.x * r.y;
        mp.expansion[XXZ] -= p_mass * r.x * r.x * r.z;
        mp.expansion[XYY] -= p_mass * r.x * r.y * r.y;
        mp.expansion[XYZ] -= p_mass * r.x * r.y * r.z;
        mp.expansion[XZZ] -= p_mass * r.x * r.z * r.z;
        mp.expansion[YYY] -= p_mass * r.y * r.y * r.y;
        mp.expansion[YYZ] -= p_mass * r.y * r.y * r.z;
        mp.expansion[YZZ] -= p_mass * r.y * r.z * r.z;
        mp.expansion[ZZZ] -= p_mass * r.z * r.z * r.z;
#endif
#if P > 3
        mp.expansion[XXXX] += p_mass * r.x * r.x * r.x * r.x;
        mp.expansion[XXXY] += p_mass * r.x * r.x * r.x * r.y;
        mp.expansion[XXXZ] += p_mass * r.x * r.x * r.x * r.z;
        mp.expansion[XXYY] += p_mass * r.x * r.x * r.y * r.y;
        mp.expansion[XXYZ] += p_mass * r.x * r.x * r.y * r.z;
        mp.expansion[XXZZ] += p_mass * r.x * r.x * r.z * r.z;
        mp.expansion[XYYY] += p_mass * r.x * r.y * r.y * r.y;
        mp.expansion[XYYZ] += p_mass * r.x * r.y * r.y * r.z;
        mp.expansion[XYZZ] += p_mass * r.x * r.y * r.z * r.z;
        mp.expansion[XZZZ] += p_mass * r.x * r.z * r.z * r.z;
        mp.expansion[YYYY] += p_mass * r.y * r.y * r.y * r.y;
        mp.expansion[YYYZ] += p_mass * r.y * r.y * r.y * r.z;
        mp.expansion[YYZZ] += p_mass * r.y * r.y * r.z * r.z;
        mp.expansion[YZZZ] += p_mass * r.y * r.z * r.z * r.z;
        mp.expansion[ZZZZ] += p_mass * r.z * r.z * r.z * r.z;
#endif
    }
    
#if P > 1
    mp.expansion[XX] *= 0.5;
    mp.expansion[YY] *= 0.5;
    mp.expansion[ZZ] *= 0.5;
#endif
#if P > 2
    mp.expansion[XXX] *= 0.1666666666666667;
    mp.expansion[YYY] *= 0.1666666666666667;
    mp.expansion[ZZZ] *= 0.1666666666666667;
    mp.expansion[XXY] *= 0.5;
    mp.expansion[XXZ] *= 0.5;
    mp.expansion[XYY] *= 0.5;
    mp.expansion[XZZ] *= 0.5;
    mp.expansion[YYZ] *= 0.5;
    mp.expansion[YZZ] *= 0.5;
#endif
#if P > 3
    mp.expansion[XXXX] *= 0.041666666666666667;
    mp.expansion[YYYY] *= 0.041666666666666667;
    mp.expansion[ZZZZ] *= 0.041666666666666667;
    mp.expansion[XXXY] *= 0.1666666666666667;
    mp.expansion[XXXZ] *= 0.1666666666666667;
    mp.expansion[XYYY] *= 0.1666666666666667;
    mp.expansion[XZZZ] *= 0.1666666666666667;
    mp.expansion[YYYZ] *= 0.1666666666666667;
    mp.expansion[YZZZ] *= 0.1666666666666667;
    mp.expansion[XXYY] *= 0.25;
    mp.expansion[XXZZ] *= 0.25;
    mp.expansion[YYZZ] *= 0.25;
    mp.expansion[XXYZ] *= 0.5;
    mp.expansion[XYYZ] *= 0.5;
    mp.expansion[XYZZ] *= 0.5;
#endif
    
//    float3 dims = mp.max - mp.min;
//    mp.size = max3(dims.x, dims.y, dims.z);
    mp.size = length(max(abs(mp.max - mp.pos), abs(mp.min - mp.pos)));

    addPowers(mp);
    
    return mp;
}
