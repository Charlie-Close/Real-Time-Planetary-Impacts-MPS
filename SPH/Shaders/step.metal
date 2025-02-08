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
    const float dt = DT_MIN1 * (*_dt);

    if (index == 0) {
        *globalTime += (*_dt);
    }
    
    float3 velocity = velocities[index];
    const float3 acceleration = accelerations[index];
    
    float h_i = h[index];
    float h1 = 1 / h_i;
    float dhdt = dhdts[index];
    h[index] = h_i * exp(h1 * dhdt * dt);
    
    if (!active[index]) {
        float density = densities[index];
        densities[index] = density * exp(-3 * h1 * dhdt * dt);
        positions[index] += dt * velocity;
        velocities[index] += dt * acceleration;
        internalEnergies[index] += dt * dInternalEnergy[index];
        return;
    }
    
    positions[index] += dt * velocity + 0.5 * acceleration * dt * dt;
    velocities[index] += dt * acceleration;
    internalEnergies[index] += dt * dInternalEnergy[index];
}
