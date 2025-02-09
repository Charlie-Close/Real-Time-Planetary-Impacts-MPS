//
//  cellsToScan.h
//  SPH
//
//  Created by Charlie Close on 15/11/2024.
//

#ifndef cellsToScan_h
#define cellsToScan_h

#include <metal_stdlib>
#include "morton.h"
using namespace metal;

struct CellToScanRange {
    int3 min;
    int3 max;
};

static inline uint cellPositionToIndex(int3 pos, int cellsPerDim) {
    int3 modded = {
        pos.x % cellsPerDim,
        pos.y % cellsPerDim,
        pos.z % cellsPerDim
    };
    
    if (modded.x < 0) {
        modded.x += cellsPerDim;
    }
    if (modded.y < 0) {
        modded.y += cellsPerDim;
    }
    if (modded.z < 0) {
        modded.z += cellsPerDim;
    }
    
    return morton3D((uint) modded.x, (uint) modded.y, (uint) modded.z);
}

static inline CellToScanRange setCellsToScanDynamic(
    float3 position,
    float cellSize,
    int cellsPerDim,
    float  h
) {
    float radius = 2.f * h;
    float3 minRange = position - float3(radius, radius, radius);
    float3 maxRange = position + float3(radius, radius, radius);
    
    float3 minNorm = minRange / cellSize;
    float3 maxNorm = maxRange / cellSize;
    
    CellToScanRange range;
    range.min = {
        (int)floor(minNorm.x),
        (int)floor(minNorm.y),
        (int)floor(minNorm.z)
    };
    range.max = {
        (int)floor(maxNorm.x),
        (int)floor(maxNorm.y),
        (int)floor(maxNorm.z),
    };

    return range;
}


#endif /* cellsToScan_h */
