//
//  density.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include "utils/kernals/kernels.h"
#include "utils/eos/eos.h"
#include "utils/cellsToScan/cellsToScan.h"
#include "utils/poles/poles.h"
#include "utils/morton.h"

constant float minRho = 1e-6;
constant float maxRho = 2.98e-2;
constant float minU = 5e+2;
constant float maxU = 8.80519e+10;

kernel void density(device float3* positions,
                    device float3* velocities,
                    device float* densities,
                    device float* internalEnergies,
                    device float* masses,
                    device float* pressures,
                    device float* h,
                    device int* materialIds,
                    device float* gradientTerms,
                    device float* speedsOfSound,
                    device float* phiGrad,
                    device float* balsara,
                    device uint2* cellData,
                    device int* cellStarts,
                    device int* cellEnds,
                    device float* cellSize,
                    device int* cellsPerDim,
                    device bool* active,
                    device int* nextActiveTime,
                    device int* globalTime,
                    device int* dt,
                    texture2d<float, access::sample> texForesite [[texture(0)]],
                    texture2d<float, access::sample> texFe [[texture(1)]],
                    uint index [[thread_position_in_grid]])
{
    if (index == 0) {
        *dt = 16 * 64;
    }
    int ind = cellData[index].y;
    if (nextActiveTime[ind] > (*globalTime)) {
        active[ind] = false;

        float u = internalEnergies[ind] * pow(10.f, 12);
        float materialId = materialIds[ind];
        float density = densities[ind];
        
        float uCoord = log(u / minU) / log(maxU / minU);
        float rhoCoord = log(density / minRho) / log(maxRho / minRho);
        float2 uv = float2(uCoord, rhoCoord);
        uv = clamp(uv, float2(0.0), float2(1.0));
        sampler textureSampler(mag_filter::linear, min_filter::linear, mip_filter::none, address::clamp_to_edge);
        float4 pc = materialId == 400 ? texForesite.sample(textureSampler, uv) : texFe.sample(textureSampler, uv);
                
        pressures[ind] = pc.x;
        speedsOfSound[ind] = pc.y;
        return;
    }
    active[ind] = true;
    
    float3 position = positions[ind];
    float h_i = h[ind];
    float h2x4 = 4 * h_i * h_i;
    float h1 = 1 / h_i;
    float mass_i = masses[ind];
    float density = 0;
    int cellsPerDimT = *cellsPerDim;
    float eps = 1;
    int count = 0;
    while (eps > 0.001 and count < 50) {
        CellToScanRange range = setCellsToScanDynamic(position, *cellSize, *cellsPerDim, h_i);
        density = 0;
        
        for (int x = range.min.x; x <= range.max.x; x++) {
            for (int y = range.min.y; y <= range.max.y; y++) {
                for (int z = range.min.z; z <= range.max.z; z++) {
                    uint cellindex = morton3D((uint) (x % cellsPerDimT),
                                              (uint) (y % cellsPerDimT),
                                              (uint) (z % cellsPerDimT));
                    uint start = cellStarts[cellindex];
                    if (start == UINT_MAX) {
                        continue;
                    }
                    uint end = cellEnds[cellindex] + 1;
                    for (uint k = start; k < end; k++) {
                        const uint j = cellData[k].y;
                        float3 x_ij = position - positions[j];
                        float r2 = length_squared(x_ij);
                        if (r2 > h2x4) {
                            continue;
                        }
                        float W_ij = W(sqrt(r2), h1);
                        float mass = masses[j];
                        density += W_ij * mass;
                    }
                }
            }
        }
                
        float newH = min(1.2348 * pow(mass_i / density, 1.f / 3), .5f);
        eps = abs(newH - h_i) / h_i;
        h_i = newH;
        h1 = 1 / h_i;
        count++;
    }
    
    
    
    CellToScanRange range = setCellsToScanDynamic(position, *cellSize, *cellsPerDim, h_i);
    density = 0;
    float density_h = 0;
    float velDiv = 0;
    float3 velCurl = { 0, 0, 0 };
    float3 v_i = velocities[ind];
    
    for (int x = range.min.x; x <= range.max.x; x++) {
        for (int y = range.min.y; y <= range.max.y; y++) {
            for (int z = range.min.z; z <= range.max.z; z++) {
                uint cellindex = morton3D((uint) (x % cellsPerDimT),
                                          (uint) (y % cellsPerDimT),
                                          (uint) (z % cellsPerDimT));

                uint start = cellStarts[cellindex];
                if (start == UINT_MAX) {
                    continue;
                }
                uint end = cellEnds[cellindex] + 1;
                for (uint k = start; k < end; k++) {
                    const uint j = cellData[k].y;
                    float3 x_ij = position - positions[j];
                    float r2 = length_squared(x_ij);
                    if (r2 > h2x4) {
                        continue;
                    }
                    float r = sqrt(r2);
                    float r1 = 1 / r;
                    float W_ij = W(r, h1);
                    float mass = masses[j];
                    density += W_ij * mass;
                    density_h -= mass * dW_dh(r, h1);
                    float3 v_j = velocities[j];
                    float3 v_ij = v_i - v_j;
                    float3 gW = gradW(x_ij, r, r1, h1);
                    velDiv += mass * dot(v_ij, gW);
                    velCurl += mass * cross(v_ij, gW);
                }
            }
        }
    }
    
    velDiv = abs(velDiv) / density;
    h[ind] = h_i;
    
    densities[ind] = density;
    gradientTerms[ind] = - 1 / (1 + (h_i / (3 * density)) * density_h);
    phiGrad[ind] = - (h_i / (3 * density)) * phiGrad[ind];
    
    float u = internalEnergies[ind] * pow(10.f, 12);
    float materialId = materialIds[ind];
    
    float uCoord = log(u / minU) / log(maxU / minU);
    float rhoCoord = log(density / minRho) / log(maxRho / minRho);
    float2 uv = float2(uCoord, rhoCoord);
    uv = clamp(uv, float2(0.0), float2(1.0));
    sampler textureSampler(mag_filter::linear, min_filter::linear, mip_filter::none, address::clamp_to_edge);
    float4 pc = materialId == 400 ? texForesite.sample(textureSampler, uv) : texFe.sample(textureSampler, uv);
            
    pressures[ind] = pc.x;
    speedsOfSound[ind] = pc.y;
    
    balsara[ind] = velDiv / (velDiv + length(velCurl) + 1e-4 * (pc.y / h_i));
}
