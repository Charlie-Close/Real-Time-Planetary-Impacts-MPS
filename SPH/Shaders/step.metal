//
//  step.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include <metal_stdlib>
using namespace metal;

kernel void step(device float3* positions,
                 device float3* velocities,
                 device float3* accelerations,
                 device float* internalEnergies,
                 device float* dInternalEnergy,
                 device bool* isAlive,
                 device int* _dt,
                 uint index [[thread_position_in_grid]])
{
    const float dt = .0001f * (*_dt);
    float3 velocity = velocities[index];
    
    if (!isAlive[index]) {
//        positions[index] += dt * velocity;
        return;
    }
    
    const float3 acceleration = accelerations[index];
    positions[index] += dt * velocity + 0.5 * acceleration * dt * dt;
    velocities[index] += dt * acceleration;
    internalEnergies[index] += dt * dInternalEnergy[index];
    
    float3 position = positions[index];
//    if (position.x > 500 or position.y > 500 or position.z > 500 or position.x < 100 or position.y < 100 or position.z < 100 or length(acceleration) > dt * 0.001) {
////    if (position.x > 380 or position.y > 380 or position.z > 380 or position.x < 240 or position.y < 240 or position.z < 240) {
//        isAlive[index] = false;
//    }
}
