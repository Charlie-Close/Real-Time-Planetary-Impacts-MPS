//
//  activate.metal
//  SPH
//
//  Created by Charlie Close on 16/03/2025.
//

#include <metal_stdlib>
using namespace metal;


kernel void activate(device bool* active,
                     device int* nextActiveTime,
                     device int& globalTime,
                     device int& _dt,
                     device bool& activateAll,
                     uint index [[thread_position_in_grid]])
{
    active[index] = activateAll || globalTime >= nextActiveTime[index];
}
