//
//  cellsToScan.h
//  SPH
//
//  Created by Charlie Close on 15/11/2024.
//

#ifndef cellsToScan_h
#define cellsToScan_h

#include <metal_stdlib>
using namespace metal;

struct CellToScanRange {
    int3 min;
    int3 max;
};

CellToScanRange setCellsToScanDynamic(
    float3 position,
    float cellSize,
    int cellsPerDim,
    float  h
);

#endif /* cellsToScan_h */
