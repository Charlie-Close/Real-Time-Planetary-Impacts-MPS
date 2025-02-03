//
//  cellsToScan.metal
//  SPH
//
//  Created by Charlie Close on 15/11/2024.
//

#include "cellsToScan.h"

CellToScanRange setCellsToScanDynamic(
    float3 position,            // Particle position
    float cellSize,            // Bounding box
    int cellsPerDim,
    float  h                    // Smoothing length
) {
    float radius = 2.f * h;
    float3 minRange = position - float3(radius, radius, radius);
    float3 maxRange = position + float3(radius, radius, radius);
    
    float3 minNorm = minRange / cellSize;
    float3 maxNorm = maxRange / cellSize;
    
    CellToScanRange range;
    range.min = {
        (int)floor(minNorm.x) % cellsPerDim,
        (int)floor(minNorm.y) % cellsPerDim,
        (int)floor(minNorm.z) % cellsPerDim
    };
    range.max = {
        (int)floor(maxNorm.x) % cellsPerDim,
        (int)floor(maxNorm.y) % cellsPerDim,
        (int)floor(maxNorm.z) % cellsPerDim,
    };
    
    if (range.min.x < 0) {
        range.min.x += cellsPerDim;
    }
    if (range.min.y < 0) {
        range.min.y += cellsPerDim;
    }
    if (range.min.z < 0) {
        range.min.z += cellsPerDim;
    }
    
    if (range.max.x < range.min.x) {
        range.max.x += cellsPerDim;
    }
    if (range.max.y < range.min.y) {
        range.max.y += cellsPerDim;
    }
    if (range.max.z < range.min.z) {
        range.max.z += cellsPerDim;
    }

    return range;
}
