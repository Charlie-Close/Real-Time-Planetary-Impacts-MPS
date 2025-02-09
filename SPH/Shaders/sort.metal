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

// Radix sort running on a GPU. This is used for finding a particle's neighbours. Has a linear time complexity.

kernel void hash(device float3* positions,
                 device float* cellSize,
                 device int* cellsPerDim,
                 device uint2* cellParticles,
                 uint index [[thread_position_in_grid]])
{
    // Get the cell index of each particle and store it in an array which we will then sort.
    float3 normalised = positions[index] / (*cellSize);
    int3 pos = {
        (int)floor(normalised.x) % (*cellsPerDim),
        (int)floor(normalised.y) % (*cellsPerDim),
        (int)floor(normalised.z) % (*cellsPerDim)
    };
    
    if (pos.x < 0) {
        pos.x += (*cellsPerDim);
    }
    if (pos.y < 0) {
        pos.y += (*cellsPerDim);
    }
    if (pos.z < 0) {
        pos.z += (*cellsPerDim);
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
    // We split particles up into blocks and produce a historgram.
    // We let the particle know where in it's blocks histogram it lies.
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

kernel void scan(device uint2* cellParticles,
                 device uint* bucketHist,
                 device uint* bucketOffset,
                 device uint* particleOffsets,
                 device uint* ittr,
                 const device int* nBlocks,
                 uint index [[thread_position_in_grid]])
{
    // Sum buckets up for each block level histogram to produce a global histogram.
    // We store this in bucket offsets (which will later actually become offsets, for now
    // it is just a histogram). We also let each block's bucket known where in the overall
    // bucket it lies.
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


kernel void sum(device uint* bucketHist,
                 device uint* bucketOffset,
                 uint index [[thread_position_in_grid]])
{
    // We now sum over all the histograms in bucket offset to produce their actuall offsets:
    // i.e. each bucket now knows where in the entire data set it starts.
    uint count = 0;
    for (int i = 0; i < SORTING_BUCKET_NUMBER; i++) {
        count += bucketOffset[i];
        bucketOffset[i] = count;
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
    // Put it all together. Now each particle knows where in the block it lies, each block knows where in
    // the bucket it lies and each bucket knows where in the array it lies, we can calculate where each particle
    // lies in the array.
    uint shift = (*ittr) * SORTING_MASK_LENGTH;
    uint bucketIndex = (cellParticles[index].x >> shift) & SORTING_BIT_MASK;
    int block = index / SORTING_BLOCK_SIZE;
    uint bOff = bucketOffset[bucketIndex]; // where the bucket lies in the array
    uint blockOffset = block == 0 ? 0 : bucketHist[(block - 1) * SORTING_BUCKET_NUMBER + bucketIndex]; // where the block lies in the bucket
    uint subBlockOffset = particleOffsets[index]; // where the particle lies in the block
    newCellParticles[bOff + subBlockOffset + blockOffset] = cellParticles[index]; // Put together.
}

kernel void initialise(device uint* cellStarts,
                       device uint* cellEnds,
                       uint index [[thread_position_in_grid]]) {
    // Now we have our sorted array, we want to know where in the array we can find each cell.
    // We store a cell start and a cell end value which we initialize as below, to show it is empty.
    cellStarts[index] = UINT_MAX;
    cellEnds[index] = 0;
}

kernel void findCellPositions(device uint2* sortedCellParticles,
                              device atomic_uint* cellStarts,
                              device atomic_uint* cellEnds,
                              uint index [[thread_position_in_grid]]) {
    int cell = sortedCellParticles[index].x;
    // Atomically min and max.
    atomic_fetch_min_explicit(&cellStarts[cell], index, memory_order_relaxed);
    atomic_fetch_max_explicit(&cellEnds[cell], index, memory_order_relaxed);
}










