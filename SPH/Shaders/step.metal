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
                 device float* drho_dts,
                 device float* h,
                 device bool* active,
                 device int* nextActiveTime,
                 device int& globalTime,
                 device int& _dt,
                 device bool* alive,
                 uint index [[thread_position_in_grid]])
{
    if (index == 0) {
        // Advance global time
        globalTime += _dt;
    }
    
    if (!alive[index]) {
        return;
    }
    
    if (active) {
        accelerations[index] = accelerations1[index];
    }
    
    // Convert dt from integer to float
    const float dt = MIN_DT * _dt;
    
    float3 position = positions[index];
    float3 velocity = velocities[index];
    float3 a = accelerations[index] + gravAccelerations[index];
    
    // Integrate forwards
    position += dt * velocity + 0.5 * a * dt * dt;
    positions[index] = position;
    velocities[index] = velocity + dt * a;
    float h_i = h[index];
    float drho_dt = drho_dts[index];
    float dhdt = - 0.3333333 * drho_dt * h_i / densities[index];
    h[index] = clamp(h_i + dhdt * dt, 1e-12, MAX_SMOOTHING_LENGTH);
    densities[index] += drho_dt * dt;
    internalEnergies[index] += dt * dInternalEnergy[index];
    
    position -= float3(BOX_CENTER);
    if (abs(position.x) > BOX_SIZE or abs(position.y) > BOX_SIZE or abs(position.z) > BOX_SIZE) {
        alive[index] = false;
    }
}
