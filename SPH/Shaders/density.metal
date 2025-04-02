//
//  density.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include "utils/kernels.h"
#include "utils/cellsToScan.h"
#include "utils/eos.h"
#include "../Parameters.h"

inline void outerProductAdd(thread float3x3& M,
                            thread const float& scale,
                            thread const float3& v,
                            thread const float3& g) {
    // unrolled outer product
    M[0][0] += scale * v.x * g.x;
    M[0][1] += scale * v.x * g.y;
    M[0][2] += scale * v.x * g.z;
    M[1][0] += scale * v.y * g.x;
    M[1][1] += scale * v.y * g.y;
    M[1][2] += scale * v.y * g.z;
    M[2][0] += scale * v.z * g.x;
    M[2][1] += scale * v.z * g.y;
    M[2][2] += scale * v.z * g.z;
}

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
                    device float3* rhoGrads,
                    device float* temperatures,
                    device uint2* cellData,
                    device int* cellStarts,
                    device int* cellEnds,
                    device uint* largeBoundCellParticles,
                    device bool* active,
                    device bool* alive,
                    device int& dt,
                    device float* pAlphaLoc,
                    device float3* accelerations,
                    device int& globalTime,
                    device int& gravNextActiveTime,
                    texture2d<float, access::sample> texIron [[texture(0)]],
                    texture2d<float, access::sample> texForesite [[texture(1)]],
                    texture2d<float, access::sample> texFe [[texture(2)]],
                    texture2d<float, access::sample> texHHe [[texture(3)]],
                    texture2d<float, access::sample> texIce [[texture(4)]],
                    texture2d<float, access::sample> texRock [[texture(5)]],
                    uint index [[thread_position_in_grid]])
{
    if (index == 0) {
        int gravDt = gravNextActiveTime - globalTime;
        if (gravDt != 0) {
            dt = gravDt;
        } else {
            dt = (int)floor(MAX_DT / MIN_DT);
        }
    }
    
    uint ind = cellData[index].y;
    
    if (!alive[ind]) {
        return;
    }
            
    // Handle case that particle is inactive.
    if (!active[ind]) {
        float u = internalEnergies[ind];
        float materialId = materialIds[ind];
        float density = densities[ind];
        float4 pc = eos(materialId, u, density, texIron, texForesite, texFe, texHHe, texIce, texRock);

        pressures[ind] = pc.x;
        speedsOfSound[ind] = pc.y;
        temperatures[ind] = pc.z;

        return;
    }
    
    float3 x_i = positions[ind];
    float mass_i = masses[ind];
    float h_i = h[ind];
    int max_h = 0;
    float hSupport = GAMMA * h_i;
    float hSupportSqrd = hSupport * hSupport;
    float h1 = 1 / h_i;
    float H1 = h1 * gamma1;
    float kernal_fac = KERNEL_CONSTANT * H1 * H1 * H1;
    float dh_kernel_fac = kernal_fac * H1;
    float density;
    float density_h;
    float omega;
    float eps = 1;
    int count = 0;
    
    float min_bound = 0;
    float max_bound = MAX_SMOOTHING_LENGTH;

    // Cache neighbours to speed up second pass
    int nNeighbours;
    uint neighbourIndex[N_NEIGHBOURS_ESTIM];
    float rs[N_NEIGHBOURS_ESTIM];

    // Find the smoothing length of our particle with the NR method.
    while (eps > DENSITY_ETA and count < MAX_DENSITY_NR_ITTERATIONS) {
        // Find which cells are within our smoothing length.
        density = 0;
        density_h = 0;
        int particlesInRange = 0;

        if (h_i < max_h and nNeighbours < N_NEIGHBOURS_ESTIM) {
            for (int i = 0; i < nNeighbours; i++) {
                const uint j = neighbourIndex[i];
                float r = rs[i];
                if (r > hSupport) {
                    // out of range, continue
                    continue;
                }
                float densityFactor = findDensityFactor(materialIds[ind], materialIds[j], temperatures[j], pressures[j], texIron, texForesite, texFe, texHHe, texIce, texRock);
                float W_ij = W_fast(r, H1, kernal_fac);
                float mass = densityFactor * masses[j];
                density += W_ij * mass;
                density_h += mass * dW_dh_fast(r, H1, kernal_fac);
                particlesInRange++;
            }
        } else {
            CellToScanRange range = setCellsToScanDynamic(x_i, h_i);
            nNeighbours = 0;
            for (int x = range.min.x; x <= range.max.x; x++) {
                for (int y = range.min.y; y <= range.max.y; y++) {
                    for (int z = range.min.z; z <= range.max.z; z++) {
                        uint cellindex = cellPositionToIndex({ x, y, z });
                        uint start = cellStarts[cellindex];
                        if (start == UINT_MAX) {
                            // Cell is empty, continue.
                            continue;
                        }
                        
                        // Check if cell is in range
                        float3 cellMin = float3(x, y, z) * CELL_WIDTH;
                        float3 cellMax = cellMin + float3(CELL_WIDTH, CELL_WIDTH, CELL_WIDTH);
                        float3 distToCell = clamp(x_i, cellMin, cellMax) - x_i; // component-wise clamp
                        float distToCell2 = dot(distToCell, distToCell);
                        if (distToCell2 > hSupportSqrd) {
                          continue;
                        }
                        
                        uint largeBoundCellIndex = cellPositionToLargeBoundIndex({ x, y, z });
                        // Loop through every particle within the cell
                        uint end = cellEnds[cellindex] + 1;
                        for (uint k = start; k < end; k++) {
                            const uint j = cellData[k].y;
                            if (largeBoundCellIndex != largeBoundCellParticles[j]) {
                                continue;
                            }
                            
                            float3 x_j = positions[j];
                            
                            float3 x_ij = x_j - x_i;
                            float r2 = length_squared(x_ij);
                            if (r2 > hSupportSqrd) {
                                // Cell out of range, continue
                                continue;
                            }
                            float r = sqrt(r2);
                            // If we have space, add to cache.
                            if (nNeighbours < N_NEIGHBOURS_ESTIM) {
                                neighbourIndex[nNeighbours] = j;
                                rs[nNeighbours] = r;
                                nNeighbours++;
                            }
                            float densityFactor = findDensityFactor(materialIds[ind], materialIds[j], temperatures[j], pressures[j], texIron, texForesite, texFe, texHHe, texIce, texRock);
                            float W_ij = W_fast(r, H1, kernal_fac);
                            float mass = densityFactor * masses[j];
                            density += W_ij * mass;
                            density_h += mass * dW_dh_fast(r, H1, dh_kernel_fac);
                            particlesInRange++;
                        }
                    }
                }
            }
        }
    
        // Estimate smoothing length based on density.
        omega = 1 + density_h * h_i / (3 * density);
        float eta_h = (RESOLUTION_ETA / h_i);
        float eta = mass_i * eta_h * eta_h * eta_h - density;
        float newH = clamp(h_i * (1 + eta / (3 * density * omega)), 1e-12, MAX_SMOOTHING_LENGTH);
        if (particlesInRange < PARTICLE_IN_RANGE_THRESHOLD) {
            newH = min(max(newH, h_i * 2), MAX_SMOOTHING_LENGTH);
        }
        
        if (newH < h_i) {
            max_bound = h_i;
        } else {
            min_bound = h_i;
        }
        
        if (newH < min_bound or newH > max_bound) {
            newH = 0.5 * (min_bound + max_bound);
        }
        
        eps = abs(newH - h_i) / h_i;
        h_i = newH;
        max_h = fmax(h_i, max_h);
        h1 = 1 / h_i;
        H1 = h1 * gamma1;
        kernal_fac = KERNEL_CONSTANT * H1 * H1 * H1;
        dh_kernel_fac = kernal_fac * H1;
        hSupport = GAMMA * h_i;
        hSupportSqrd = hSupport * hSupport;
        count++;
    }

    // Now we have our smoothing length. We can go and calculate this stuff:
    float velDiv = 0;
    float3 velCurl(0);
    float3 v_i = velocities[ind];
    float3 rhoGrad(0);
    float3x3 M(0);
    float accDiv = 0;
    float unscaledDensity = 0;
    
    // If we haven't overflowed, we can use our cache to speed things up.
    if (nNeighbours < N_NEIGHBOURS_ESTIM) {
        for (int i = 0; i < nNeighbours; i++) {
            const uint j = neighbourIndex[i];
            float r = rs[i];
            if (r > hSupport) {
                // Particle is out of range - skip it.
                continue;
            }
            // Apply our equations.
            float mass = masses[j];
            if (r > 1e-12) {
                float r1 = 1 / r;
                float3 v_j = velocities[j];
                float3 v_ij = v_j - v_i;
                float3 x_ij = positions[j] - x_i;
                float3 gW = gradW_fast(x_ij, r, r1, H1, dh_kernel_fac);
                velDiv += mass * dot(v_ij, gW);
                velCurl += mass * cross(v_ij, gW);
                outerProductAdd(M, -mass, v_ij, gW);
                accDiv += mass * dot(accelerations[j] - accelerations[ind], gW);
                
                // Non physical, just for rendering purposes
                rhoGrad += 0.25 * mass * (density * gW + 3 * rhoGrads[j] * W_fast(r, H1, kernal_fac));
            }
            unscaledDensity += mass * W_fast(r, H1, kernal_fac);
        }
    } else {
        // Our cache has overflowed. Regrettably we need to loop though surrounding cells.
        CellToScanRange range = setCellsToScanDynamic(x_i, h_i);
        for (int x = range.min.x; x <= range.max.x; x++) {
            for (int y = range.min.y; y <= range.max.y; y++) {
                for (int z = range.min.z; z <= range.max.z; z++) {
                    uint cellindex = cellPositionToIndex({ x, y, z });

                    uint start = cellStarts[cellindex];
                    if (start == UINT_MAX) {
                        continue;
                    }
                    
                    // Check if cell is in range
                    float3 cellMin = float3(x, y, z) * CELL_WIDTH;
                    float3 cellMax = cellMin + float3(CELL_WIDTH, CELL_WIDTH, CELL_WIDTH);
                    float3 distToCell = clamp(x_i, cellMin, cellMax) - x_i; // component-wise clamp
                    float distToCell2 = dot(distToCell, distToCell);
                    if (distToCell2 > hSupportSqrd) {
                      continue;
                    }
                    
                    uint end = cellEnds[cellindex] + 1;
                    for (uint k = start; k < end; k++) {
                        const uint j = cellData[k].y;
                        if (!alive[j]) {
                            continue;
                        }
                        
                        float3 x_j = positions[j];
                        
                        float3 x_ij = x_j - x_i;
                        float r2 = length_squared(x_ij);
                        if (r2 > hSupport) {
                            // Particle is out of range - skip it.
                            continue;
                        }
                        // Apply our equations.
                        float r = sqrt(r2);
                        float mass = masses[j];
                        if (r > 1e-12) {
                            float r1 = 1 / r;
                            float3 v_j = velocities[j];
                            float3 v_ij = v_j - v_i;
                            float3 gW = gradW_fast(x_ij, r, r1, H1, dh_kernel_fac);
                            velDiv += mass * dot(v_ij, gW);
                            velCurl += mass * cross(v_ij, gW);
                            outerProductAdd(M, - mass, v_ij, gW);
                            accDiv += mass * dot(accelerations[j] - accelerations[ind], gW);
                            
                            // Non physical, just for rendering purposes
                            rhoGrad += 0.25 * mass * (density * gW + 3 * rhoGrads[j] * W_fast(r, H1, kernal_fac));
                        }
                        unscaledDensity += mass * W_fast(r, H1, kernal_fac);
                    }
                }
            }
        }
    }
            
    float sum = 0.0f;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            sum += M.columns[i][j] * M.columns[j][i];
        }
    }

    rhoGrads[ind] = rhoGrad / unscaledDensity;
    velDiv /= unscaledDensity;
    float absCurl = length(velCurl) / unscaledDensity;

    // Update our smoothing length, density and gradient terms
    h[ind] = h_i;
    densities[ind] = density;

    // Apply our equations of state to update our pressure and sound speed
    float u = internalEnergies[ind];
    int materialId = materialIds[ind];

    float4 pc = eos(materialId, u, density, texIron, texForesite, texFe, texHHe, texIce, texRock);

    pressures[ind] = pc.x;
    speedsOfSound[ind] = pc.y;
    temperatures[ind] = pc.z;
    
    // Calculate our balsara switch.
    balsara[ind] = abs(velDiv) / (abs(velDiv) + absCurl + 1e-4 * (pc.y / h_i));
    gradientTerms[ind] = 1 / omega;
    
    if (nNeighbours > 5) {
        float S = 10 * h_i * h_i * max(0.f, -1.f * (accDiv - sum) / unscaledDensity);
        pAlphaLoc[ind] = clamp(ALPHA_MAX * S / (pc.y * pc.y + S), ALPHA_MIN, ALPHA_MAX);
    } else {
        pAlphaLoc[ind] = ALPHA_MIN;
    }
}
