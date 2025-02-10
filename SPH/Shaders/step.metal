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
                 device float* densities,
                 device float* internalEnergies,
                 device float* dInternalEnergy,
                 device float* dhdts,
                 device float* h,
                 device bool* active,
                 device int* nextActiveTime,
                 device int* globalTime,
                 device int* _dt,
                 uint index [[thread_position_in_grid]])
{
    // Convert dt from integer to float
    const float dt = DT_MIN1 * (*_dt);

    if (index == 0) {
        // Advance global time
        *globalTime += (*_dt);
    }
    
    float3 velocity = velocities[index];
    const float3 acceleration = accelerations[index];
    
    // Integrate forwards
    float h_i = h[index];
    float h1 = 1 / h_i;
    float dhdt = dhdts[index];
    h[index] = min(h_i * exp(h1 * dhdt * dt), MAX_SMOOTHING_LENGTH);
    densities[index] *= exp(-3 * h1 * dhdt * dt);
    velocities[index] = velocity + dt * acceleration;
    internalEnergies[index] += dt * dInternalEnergy[index];
    positions[index] += dt * velocity + 0.5 * acceleration * dt * dt;
}
