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
                 device bool* alive,
                 uint index [[thread_position_in_grid]])
{
    if (index == 0) {
        // Advance global time
        *globalTime += (*_dt);
    }
    
    // Convert dt from integer to float
    const float dt = DT_MIN1 * (*_dt);
    
    if (!alive[index]) {
        return;
    }
    
    float3 position = positions[index];
    float3 velocity = velocities[index];
    const float3 acceleration = accelerations[index];
    
    // Integrate forwards
    position += dt * velocity + 0.5 * acceleration * dt * dt;
    positions[index] = position;
    velocities[index] = velocity + dt * acceleration;
    float h_i = h[index];
    float h1 = 1 / h_i;
    float dhdt = dhdts[index];
    h[index] = clamp(h_i * exp(h1 * dhdt * dt), 1e-12, MAX_SMOOTHING_LENGTH);
    densities[index] *= exp(-3 * h1 * dhdt * dt);
    internalEnergies[index] += dt * dInternalEnergy[index];
    
    position -= float3(BOX_CENTER);
    if (abs(position.x) > BOX_SIZE or abs(position.y) > BOX_SIZE or abs(position.z) > BOX_SIZE) {
        alive[index] = false;
    }
}
