//
//  sort.metal
//  SPH
//
//  Created by Charlie Close on 31/01/2025.
//

#include <metal_stdlib>
#include "utils/morton.h"
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
    int blockSize = 256;
    uint shift = (*ittr) * 8;
    uint localHist[256];
    for (int i = 0; i < 256; i++) {
        localHist[i] = 0;
    }
    
    for (int i = blockSize * index; i < min(blockSize * ((int)index + 1), (*nParticles)); i++) {
        uint bucketIndex = (cellParticles[i].x >> shift) & 0xFF;
        particleOffsets[i] = localHist[bucketIndex];
        localHist[bucketIndex]++;
    }
    
    for (int i = 0; i < 256; i++) {
        bucketHist[index * 256 + i] = localHist[i];
    }
}


kernel void sum(device uint* bucketHist,
                 device uint* bucketOffset,
                 uint index [[thread_position_in_grid]])
{
    uint count = 0;
    for (int i = 0; i < 256; i++) {
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
        count += bucketHist[i * 256 + index];
        bucketHist[i * 256 + index] = count;
    }
//    bucketOffset[index] = count;

    if (index < 255) {
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
    uint shift = (*ittr) * 8;
    uint bucketIndex = (cellParticles[index].x >> shift) & 0xFF;
    int block = index / 256;
    uint bOff = bucketOffset[bucketIndex];
    uint blOffset = block == 0 ? 0 : bucketHist[(block - 1) * 256 + bucketIndex];
    uint pOff = particleOffsets[index];
    newCellParticles[bOff + pOff + blOffset] = cellParticles[index];
}



//kernel void scan(device uint2* cellParticles,
//                 device uint* bucketHist,
//                 device uint* bucketOffset,
//                 device uint* particleOffsets,
//                 device uint* ittr,
//                 const device int* nParticles,
//                 uint index [[thread_position_in_grid]])
//{
//    int blockSize = 256;
//    uint shift = (*ittr) * 8;
//    int sumHist[256];
//    for (int i = 0; i < 256; i++) {
//        sumHist[i] = 0;
//    }
//    for (uint i = 0; i < index; i++) {
//        for (int j = 0; j < 256; j++) {
//            sumHist[j] += bucketHist[i * 256 + j];
//        }
//    }
//    
//    for (int i = blockSize * index; i < min(blockSize * ((int)index + 1), (*nParticles)); i++) {
//        uint bucketIndex = (cellParticles[i].x >> shift) & 0xFF;
//        particleOffsets[i] += sumHist[bucketIndex];
//    }
//    if (blockSize * ((int)index + 1) >= (*nParticles)) {
//        int count = 0;
//        for (int i = 0; i < 256; i++) {
//            bucketOffset[i] = count;
//            count += sumHist[i] + bucketHist[index * 256 + i];
//        }
//    }
//}
//
//kernel void sort(device uint2* cellParticles,
//                 device uint2* newCellParticles,
//                 device uint* bucketOffset,
//                 device uint* particleOffsets,
//                 device uint* ittr,
//                 uint index [[thread_position_in_grid]])
//{
//    uint shift = (*ittr) * 8;
//    uint bucketIndex = (cellParticles[index].x >> shift) & 0xFF;
//    uint bOff = bucketOffset[bucketIndex];
//    uint pOff = particleOffsets[index];
//    newCellParticles[bOff + pOff] = cellParticles[index];
//}









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










