//
//  density.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include "utils/kernels.h"
#include "utils/cellsToScan.h"
#include "utils/morton.h"
#include "../Parameters.h"

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
                    device float* balsara,
                    device uint2* cellData,
                    device int* cellStarts,
                    device int* cellEnds,
                    device float* cellSize,
                    device int* cellsPerDim,
                    device bool* active,
                    device int* dt,
                    texture2d<float, access::sample> texForesite [[texture(0)]],
                    texture2d<float, access::sample> texFe [[texture(1)]],
                    uint index [[thread_position_in_grid]])
{
    if (index == 0) {
        *dt = DT0 * DT_MIN;
    }
    int ind = cellData[index].y;
    
    // Handle case that particle is inactive.
    if (!active[ind]) {
        float u = internalEnergies[ind];
        float materialId = materialIds[ind];
        float density = densities[ind];
        
        float uCoord = log(u / ANEOS_MIN_U) / log(ANEOS_MAX_U / ANEOS_MIN_U);
        float rhoCoord = log(density / ANEOS_MIN_RHO) / log(ANEOS_MAX_RHO / ANEOS_MIN_RHO);
        float2 uv = float2(uCoord, rhoCoord);
        uv = clamp(uv, float2(0.0), float2(1.0));
        sampler textureSampler(mag_filter::linear, min_filter::linear, mip_filter::none, address::clamp_to_edge);
        float4 pc = materialId == 400 ? texForesite.sample(textureSampler, uv) : texFe.sample(textureSampler, uv);
                
        pressures[ind] = pc.x;
        speedsOfSound[ind] = pc.y;
        return;
    }
    
    float3 position = positions[ind];
    float h_i = h[ind];
    float h2x4 = 4 * h_i * h_i;
    float h1 = 1 / h_i;
    float mass_i = masses[ind];
    float density = 0;
    int cellsPerDimT = *cellsPerDim;
    float eps = 1;
    int count = 0;

    // Cache neighbours to speed up second pass
    int nNeighbours;
    uint neighbourIndex[N_NEIGHBOURS_ESTIM];
        
    // Find the smoothing length of our particle with the NR method.
    while (eps > DENSITY_ETA and count < MAX_DENSITY_NR_ITTERATIONS) {
        // Find which cells are within our smoothing length.
        CellToScanRange range = setCellsToScanDynamic(position, *cellSize, *cellsPerDim, h_i);
        density = 0;
        nNeighbours = 0;
        
        for (int x = range.min.x; x <= range.max.x; x++) {
            for (int y = range.min.y; y <= range.max.y; y++) {
                for (int z = range.min.z; z <= range.max.z; z++) {
                    uint cellindex = morton3D((uint) (x % cellsPerDimT),
                                              (uint) (y % cellsPerDimT),
                                              (uint) (z % cellsPerDimT));
                    
                    uint start = cellStarts[cellindex];
                    if (start == UINT_MAX) {
                        // Cell is empty, continue.
                        continue;
                    }
                    
                    // Loop through every particle within the cell
                    uint end = cellEnds[cellindex] + 1;
                    for (uint k = start; k < end; k++) {
                        const uint j = cellData[k].y;
                        float3 x_ij = position - positions[j];
                        float r2 = length_squared(x_ij);
                        if (r2 > h2x4) {
                            // Cell out of range, continue
                            continue;
                        }
                        // If we have space, add to cache.
                        if (nNeighbours < N_NEIGHBOURS_ESTIM) {
                            neighbourIndex[nNeighbours] = j;
                            nNeighbours++;
                        }
                        float W_ij = W(sqrt(r2), h1);
                        float mass = masses[j];
                        density += W_ij * mass;
                    }
                    
                }
            }
        }
         
        // Estimate smoothing length based on density.
        float newH = min(1.2348 * pow(mass_i / density, 1.f / 3), MAX_SMOOTHING_LENGTH);
        // eps checks if we are converging.
        eps = abs(newH - h_i) / h_i;
        h_i = newH;
        h1 = 1 / h_i;
        h2x4 = 4 * h_i * h_i;
        count++;
    }
    
    // Now we have our smoothing length. We can go and calculate this stuff:
    density = 0;
    float density_h = 0;
    float velDiv = 0;
    float3 velCurl = { 0, 0, 0 };
    float3 v_i = velocities[ind];
    // If we haven't overflowed, we can use our cache to speed things up.
    if (nNeighbours < N_NEIGHBOURS_ESTIM) {
        for (int i = 0; i < nNeighbours; i++) {
            int j = neighbourIndex[i];
            float3 x_ij = position - positions[j];
            float r2 = length_squared(x_ij);
            if (r2 > h2x4) {
                // Particle is out of range - skip it.
                continue;
            }
            // Apply our equations.
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
    } else {
        // Our cache has overflowed. Regrettably we need to loop though surrounding cells.
        CellToScanRange range = setCellsToScanDynamic(position, *cellSize, *cellsPerDim, h_i);
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
                            // Particle is out of range - skip it.
                            continue;
                        }
                        // Apply our equations.
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
    }
    
    // Normalise the divergence and curl for density
    velDiv = abs(velDiv) / density;
    float absCurl = length(velCurl) / density;
    
    // Update our smoothing length, density and gradient terms
    h[ind] = h_i;
    densities[ind] = density;
    gradientTerms[ind] = - 1 / (1 + (h_i / (3 * density)) * density_h);
    
    // Apply our equations of state to update our pressure and sound speed
    float u = internalEnergies[ind];
    float materialId = materialIds[ind];
    
    float uCoord = log(u / ANEOS_MIN_U) / log(ANEOS_MAX_U / ANEOS_MIN_U);
    float rhoCoord = log(density / ANEOS_MIN_RHO) / log(ANEOS_MAX_RHO / ANEOS_MIN_RHO);
    float2 uv = float2(uCoord, rhoCoord);
    uv = clamp(uv, float2(0.0), float2(1.0));
    sampler textureSampler(mag_filter::linear, min_filter::linear, mip_filter::none, address::clamp_to_edge);
    float4 pc = materialId == 400 ? texForesite.sample(textureSampler, uv) : texFe.sample(textureSampler, uv);
    pressures[ind] = pc.x;
    speedsOfSound[ind] = pc.y;
    
    // Calculate our balsara switch.
    balsara[ind] = velDiv / (velDiv + absCurl + 1e-4 * (pc.y / h_i));
}
