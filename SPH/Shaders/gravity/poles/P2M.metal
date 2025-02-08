//
//  P2M.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include "poles.h"
#include <metal_stdlib>
using namespace metal;

Multipole P2M(device int* treeStructure, device float* masses, device float3* positions, device float* grav, device bool* active, device int* nextActiveTime, device int* globalTime, int treePointer) {
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
    
    for (int i = start; i < end; i++) {
        int p = treeStructure[i];
        float p_mass = masses[p];
        float3 p_position = positions[p];
        int isActive = (*globalTime) <= nextActiveTime[p];
        if (!mp.active and isActive) {
            mp.active = true;
        }
        mp.max = max(mp.max, p_position);
        mp.min = min(mp.min, p_position);
        
        mp.expansion[M] += p_mass;
        mp.pos += p_position * p_mass;
        if (isActive) {
            mp.minGrav = min(mp.minGrav, grav[p]);
        }
    }
    
    for (int i = start; i < end; i++) {
        int p = treeStructure[i];
        active[p] = mp.active;
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
    
    float3 dims = mp.max - mp.min;
    mp.size = max3(dims.x, dims.y, dims.z);
    
    addPowers(mp);
    
    return mp;
}
