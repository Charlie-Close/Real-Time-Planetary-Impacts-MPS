//
//  sort.metal
//  SPH
//
//  Created by Charlie Close on 31/01/2025.
//

#include <metal_stdlib>
#include "utils/morton.h"
#include "../Parameters.h"
using namespace metal;

kernel void hash(device float3* positions,
                 device float* cellSize,
                 device int* cellsPerDim,
                 device uint2* cellParticles,
                 uint index [[thread_position_in_grid]])
{
    float3 normalised = positions[index] / (*cellSize);
    int3 pos = {
        (int)floor(normalised.x) % (*cellsPerDim),
        (int)floor(normalised.y) % (*cellsPerDim),
        (int)floor(normalised.z) % (*cellsPerDim)
    };
    
    if (normalised.x < 0) {
        normalised.x += (*cellsPerDim);
    }
    if (normalised.y < 0) {
        normalised.y += (*cellsPerDim);
    }
    if (normalised.z < 0) {
        normalised.z += (*cellsPerDim);
    }
    // Interleave the bits to get a Morton index:
    uint cellIndex = morton3D((uint) pos.x,
                              (uint) pos.y,
                              (uint) pos.z);
    cellParticles[index] = { cellIndex, index };
}

kernel void hist(device uint2* cellParticles,
                 device uint* bucketHist,
                 device uint* particleOffsets,
                 device uint* ittr,
                 const device int* nParticles,
                 uint index [[thread_position_in_grid]])
{
    uint shift = (*ittr) * SORTING_MASK_LENGTH;
    uint localHist[SORTING_BUCKET_NUMBER];
    for (int i = 0; i < SORTING_BUCKET_NUMBER; i++) {
        localHist[i] = 0;
    }
    
    for (int i = SORTING_BLOCK_SIZE * index; i < min(SORTING_BLOCK_SIZE * ((int)index + 1), (*nParticles)); i++) {
        uint bucketIndex = (cellParticles[i].x >> shift) & SORTING_BIT_MASK;
        particleOffsets[i] = localHist[bucketIndex];
        localHist[bucketIndex]++;
    }
    
    for (int i = 0; i < SORTING_BUCKET_NUMBER; i++) {
        bucketHist[index * SORTING_BUCKET_NUMBER + i] = localHist[i];
    }
}


kernel void sum(device uint* bucketHist,
                 device uint* bucketOffset,
                 uint index [[thread_position_in_grid]])
{
    uint count = 0;
    for (int i = 0; i < SORTING_BUCKET_NUMBER; i++) {
        count += bucketOffset[i];
        bucketOffset[i] = count;
    }
}

kernel void scan(device uint2* cellParticles,
                 device uint* bucketHist,
                 device uint* bucketOffset,
                 device uint* particleOffsets,
                 device uint* ittr,
                 const device int* nBlocks,
                 uint index [[thread_position_in_grid]])
{
    uint count = 0;
    for (int i = 0; i < (*nBlocks); i++) {
        count += bucketHist[i * SORTING_BUCKET_NUMBER + index];
        bucketHist[i * SORTING_BUCKET_NUMBER + index] = count;
    }

    if (index < SORTING_BIT_MASK) {
        bucketOffset[index + 1] = count;
    } else {
        bucketOffset[0] = 0;
    }
}

kernel void sort(device uint2* cellParticles,
                 device uint2* newCellParticles,
                 device uint* bucketHist,
                 device uint* bucketOffset,
                 device uint* particleOffsets,
                 device uint* ittr,
                 uint index [[thread_position_in_grid]])
{
    uint shift = (*ittr) * SORTING_MASK_LENGTH;
    uint bucketIndex = (cellParticles[index].x >> shift) & SORTING_BIT_MASK;
    int block = index / SORTING_BLOCK_SIZE;
    uint bOff = bucketOffset[bucketIndex];
    uint blockOffset = block == 0 ? 0 : bucketHist[(block - 1) * SORTING_BUCKET_NUMBER + bucketIndex];
    uint subBlockOffset = particleOffsets[index];
    newCellParticles[bOff + subBlockOffset + blockOffset] = cellParticles[index];
}

kernel void initialise(device uint* cellStarts,
                       device uint* cellEnds,
                       uint index [[thread_position_in_grid]]) {
    cellStarts[index] = UINT_MAX;
    cellEnds[index] = 0;
}

kernel void findCellPositions(device uint2* sortedCellParticles,
                              device atomic_uint* cellStarts,
                              device atomic_uint* cellEnds,
                              uint index [[thread_position_in_grid]]) {
    int cell = sortedCellParticles[index].x;
    atomic_fetch_min_explicit(&cellStarts[cell], index, memory_order_relaxed);
    atomic_fetch_max_explicit(&cellEnds[cell], index, memory_order_relaxed);
}










