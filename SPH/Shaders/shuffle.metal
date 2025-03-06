//
//  shuffle.metal
//  SPH
//
//  Created by Charlie Close on 11/02/2025.
//

#include <metal_stdlib>
using namespace metal;


kernel void shuffleFloat(device float* unsorted,
                         device float* sorted,
                         device uint2* cellParticles,
                         uint index [[thread_position_in_grid]])
{
    sorted[index] = unsorted[cellParticles[index].y];
}

kernel void shuffleFloat3(device float3* unsorted,
                          device float3* sorted,
                          device uint2* cellParticles,
                          uint index [[thread_position_in_grid]])
{
    sorted[index] = unsorted[cellParticles[index].y];
}

kernel void shuffleInt(device int* unsorted,
                       device int* sorted,
                       device uint2* cellParticles,
                       uint index [[thread_position_in_grid]])
{
    sorted[index] = unsorted[cellParticles[index].y];
}

kernel void shuffleBool(device bool* unsorted,
                       device bool* sorted,
                       device uint2* cellParticles,
                       uint index [[thread_position_in_grid]])
{
    sorted[index] = unsorted[cellParticles[index].y];
}

kernel void inverseArray(device uint2* cellParticles,
                         device uint2* cellParticlesInverse,
                         uint index [[thread_position_in_grid]]) {
    cellParticlesInverse[cellParticles[index].y].y = index;
}

kernel void shuffleTree(device int* treeStructure,
                        device uint2* cellParticlesInverse,
                        device int* pointers,
                        uint index [[thread_position_in_grid]])
{
    int treePointer = pointers[index];
    int nParticles = treeStructure[treePointer];
    
    if (nParticles != 0) {
        // We are looking at a leaf so shuffle the data.
        int start = treePointer + 2;
        int end = start + nParticles;
        for (int i = start; i < end; i++) {
            int particlePointer = treeStructure[i];
            treeStructure[i] = cellParticlesInverse[particlePointer].y;
        }
    }
}

