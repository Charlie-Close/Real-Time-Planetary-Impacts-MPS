//
//  step.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include <metal_stdlib>
#include "../Parameters.h"
using namespace metal;

kernel void step(device float3* positions,
                 device float3* velocities,
                 device float3* accelerations,
                 device float3* accelerations1,
                 device float3* gravAccelerations,
                 device float* densities,
                 device float* internalEnergies,
                 device float* dInternalEnergy,
                 device float* dh_dts,
                 device float* h,
                 device bool* active,
                 device int* nextActiveTime,
                 device int& globalTime,
                 device int& _dt,
                 device bool* alive,
                 device float* alpha,
                 device float* da_dt,
                 device float* alphaLoc,
                 uint index [[thread_position_in_grid]])
{
    if (index == 0) {
        // Advance global time
        globalTime += _dt;
    }
    
    if (!alive[index]) {
        return;
    }
    
    if (active[index]) {
        accelerations[index] = accelerations1[index];
    }

    // Convert dt from integer to float
    const float dt = MIN_DT * _dt;
    
    float3 a = accelerations[index] + gravAccelerations[index];
    
    // Integrate forwards
    positions[index] += dt * velocities[index] + 0.5 * a * dt * dt;
    velocities[index] += dt * a;
    float h_i = h[index];
    h[index] = h_i * exp((1 / h_i) * dh_dts[index] * dt);
    densities[index] *= exp(-(3.f / h_i) * dh_dts[index] * dt);
    dh_dts[index] *= h[index] / h_i;
    internalEnergies[index] = max(internalEnergies[index] + dt * dInternalEnergy[index], 0.f);
    float alpha_i = alpha[index];
    float alpha_loc = alphaLoc[index];
    float dadt = da_dt[index];
    float tau = (alpha_loc - alpha_i) / dadt;
    if (alpha_i != alpha_loc) {
        alpha_i = alpha_loc + (alpha_i - alpha_loc) * exp(dadt);
        da_dt[index] = (alpha_loc - alpha_i) / tau;
        alpha[index] = clamp(alpha_i, ALPHA_MIN, ALPHA_MAX);
    }
    
    float3 translated_pos = positions[index] - float3(BOX_CENTER);
    if (abs(translated_pos.x) > BOX_SIZE or abs(translated_pos.y) > BOX_SIZE or abs(translated_pos.z) > BOX_SIZE) {
        alive[index] = false;
    }
}
