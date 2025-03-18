//
//  accelerations.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include "utils/kernels.h"
#include "utils/cellsToScan.h"
#include "utils/eos.h"
#include "../Parameters.h"

kernel void acceleration(device float3* positions,
                         device float3* velocities,
                         device float3* accelerations,
                         device float3* accelerations1,
                         device float3* grav_accelerations,
                         device float* densities,
                         device float* internalEnergies,
                         device float* masses,
                         device float* h,
                         device float* gradientTerms,
                         device float* speedsOfSound,
                         device float* dInternalEnergy,
                         device float* balsara,
                         device float* drho_dts,
                         device uint2* cellData,
                         device int* cellStarts,
                         device int* cellEnds,
                         device uint* largeBoundCellParticles,
                         device bool* active,
                         device bool* alive,
                         device int* nextActiveTime,
                         device int& globalTime,
                         device atomic_int* dt,
                         device float* pressures,
                         device int* materialIds,
                         texture2d<float, access::sample> texIron [[texture(0)]],
                         texture2d<float, access::sample> texForesite [[texture(1)]],
                         texture2d<float, access::sample> texFe [[texture(2)]],
                         texture2d<float, access::sample> texHHe [[texture(3)]],
                         texture2d<float, access::sample> texIce [[texture(4)]],
                         texture2d<float, access::sample> texRock [[texture(5)]],
                         uint index [[thread_position_in_grid]],
                         uint localId [[thread_position_in_threadgroup]],
                         uint threadsPerGroup [[threads_per_threadgroup]])
{
    uint ind = cellData[index].y;
    
    threadgroup uint tgIndexes[THREADGROUP_CACHE_SIZE];
    threadgroup packed_float3 tgXs[THREADGROUP_CACHE_SIZE];
    
    for (uint i = localId; i < THREADGROUP_CACHE_SIZE; i += threadsPerGroup) {
        tgIndexes[i] = UINT_MAX;
    }
    
    threadgroup_barrier(mem_flags::mem_threadgroup);
    
    if (!alive[ind]) {
        return;
    }
    
    int key = ind & (THREADGROUP_CACHE_SIZE - 1);
    tgIndexes[key] = ind;
    
    threadgroup_barrier(mem_flags::mem_threadgroup);
    float3 x_i = positions[ind];
    
    if (tgIndexes[key] == ind) {
        tgXs[key] = x_i;
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);
    
    // If we are inactive this timestep, we skip
    if (!active[ind]) {
        int integerDt = max(nextActiveTime[ind] - globalTime, 1);
        atomic_fetch_min_explicit(&(*dt), integerDt, memory_order_relaxed);
        return;
    }
    
    // Get all of our relevant terms
    float rho_i = densities[ind];
    float fact_i = pressures[ind] * gradientTerms[ind] / (rho_i * rho_i);
    float3 v_i = velocities[ind];
    float c_i = speedsOfSound[ind];
    float h_i = h[ind];
    float B_i = balsara[ind];
    float hi1 = 1 / h_i;
    float hi2x4 = h_i * h_i * 4;
    

    // Accumulators
    float3 dv_dt(0);
    float du_dt = 0;
    float drho_dt = 0;
    
    // Get our surrounding cells
    CellToScanRange range = setCellsToScanDynamic(x_i, h_i);
    float v_sigi = 2 * c_i;
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
                if (distToCell2 > hi2x4) {
                  continue;
                }
                
                
                uint largeBoundCellIndex = cellPositionToLargeBoundIndex({ x, y, z });
                uint end = cellEnds[cellindex] + 1;
                // For each particle within the cell
                for (uint k = start; k < end; k++) {
                    const uint j = cellData[k].y;
                    if (largeBoundCellIndex != largeBoundCellParticles[j]) {
                        continue;
                    }

                    // Check threadgroup cache
                    uint key = j & (THREADGROUP_CACHE_SIZE - 1);
                    uint pInd = tgIndexes[key];
                    float3 x_j = pInd == j ? tgXs[key] : positions[j];
                    
                    float3 x_ij = x_i - x_j;
                    float r = length(x_ij);
                    float h_j = h[j];
                    if ((r > h_i + h_i and r > h_j + h_j) or r < 1e-12) {
                        // Cell is out of range or on top of us, skip.
                        continue;
                    }
                    // Fetch and calculate info about particle
                    float r1 = 1 / r;
                    float m_j = masses[j];
                    float rho_j = densities[j];
                    float fact_j = pressures[j] * gradientTerms[j] / (rho_j * rho_j);
                    float3 v_j = velocities[j];
                    float c_j = speedsOfSound[j];
                    float hj1 = 1 / h_j;
                    float3 v_ij = v_i - v_j;
                    float rho_ij = 0.5 * (rho_i + rho_j);
                    const float v_dot_r = dot(v_ij, x_ij);
                    float mu_ij = v_dot_r < 0 ? v_dot_r * r1 : 0;
                    float v_sigij = c_i + c_j - VISCOSITY_BETA * mu_ij;
                    v_sigi = max(v_sigij, v_sigi);
                    
                    // Viscosity term
                    float B_j = balsara[j];
                    float B_ij = 0.5 * (B_i + B_j);
                    float nu_ij = - 0.5 * VISCOSITY_ALPHA * B_ij * mu_ij * v_sigij / rho_ij;
                    
                    float3 gW_i = gradW(x_ij, r, r1, hi1);
                    float3 gW_j = gradW(x_ij, r, r1, hj1);
                    float3 gW_ij = 0.5 * (gW_i + gW_j);
                    
                    // Acceleration terms
                    float3 t_i = fact_i * gW_i;
                    float3 t_j = fact_j * gW_j;
                    float3 t_visc = nu_ij * gW_ij;
                    
                    // Accumulate
                    dv_dt -= m_j * (t_i + t_j + t_visc);
                    du_dt += m_j * (fact_i * dot(v_ij, gW_i) + 0.5 * nu_ij * dot(v_ij, gW_ij));
                    drho_dt += (m_j / rho_j) * dot(v_ij, gW_i);
                }
            }
        }
    }

    // Time stepping criterion.
    float dtCFL = 2.0f * CFL * GAMMA * h_i / v_sigi;
    float goaldt = clamp(dtCFL, MIN_DT, MAX_DT);
    float roundedDt = MIN_DT;
    while (roundedDt < goaldt and roundedDt < MAX_DT) {
        roundedDt += roundedDt;
    }
    roundedDt = 0.5 * roundedDt;
//    roundedDt = min(roundedDt, 0.5 * MIN_DT * sqrt(globalTime + 1.f));
    int integerDt = (int)max(floor(roundedDt / MIN_DT), 1.f);
    
    nextActiveTime[ind] = globalTime + integerDt;
    // dt is an integer as atomic min doesn't work with floats.
    atomic_fetch_min_explicit(&(*dt), integerDt, memory_order_relaxed);
    dInternalEnergy[ind] = du_dt;
    drho_dt *= rho_i;
    drho_dts[ind] = drho_dt;
    accelerations[ind] = dv_dt;
}

kernel void accelerationStep(device float3* positions,
                             device float3* velocities,
                             device float3* accelerations,
                             device float3* accelerations1,
                             device float3* grav_accelerations,
                             device float* densities,
                             device float* internalEnergies,
                             device float* masses,
                             device float* h,
                             device float* gradientTerms,
                             device float* speedsOfSound,
                             device float* dInternalEnergy,
                             device float* balsara,
                             device float* drho_dts,
                             device uint2* cellData,
                             device int* cellStarts,
                             device int* cellEnds,
                             device uint* largeBoundCellParticles,
                             device bool* active,
                             device bool* alive,
                             device int* nextActiveTime,
                             device int& globalTime,
                             device atomic_int* dt,
                             device float* pressures,
                             device int* materialIds,
                             texture2d<float, access::sample> texIron [[texture(0)]],
                             texture2d<float, access::sample> texForesite [[texture(1)]],
                             texture2d<float, access::sample> texFe [[texture(2)]],
                             texture2d<float, access::sample> texHHe [[texture(3)]],
                             texture2d<float, access::sample> texIce [[texture(4)]],
                             texture2d<float, access::sample> texRock [[texture(5)]],
                             uint index [[thread_position_in_grid]],
                             uint localId [[thread_position_in_threadgroup]],
                             uint threadsPerGroup [[threads_per_threadgroup]])
{
    uint ind = cellData[index].y;
    
    // If we are inactive this timestep, we skip
    if (!active[ind]) {
        return;
    }
    
    float ownDt = MIN_DT * (nextActiveTime[ind] - globalTime);
        
    // Get all of our relevant terms
    float3 a0 = accelerations[ind] + grav_accelerations[ind];
    float3 x_i = positions[ind] + ownDt * velocities[ind] + 0.5 * ownDt * ownDt * a0;
    float rho_i = densities[ind] + ownDt * drho_dts[ind];
    float3 v_i = velocities[ind] + ownDt * a0;
    float u = internalEnergies[ind] + ownDt * dInternalEnergy[ind];
    float4 pc = eos(materialIds[ind], u, rho_i, texIron, texForesite, texFe, texHHe, texIce, texRock);
    float fact_i = gradientTerms[ind] * pc.x / (rho_i * rho_i);
    float c_i = pc.y;
    float h_i = h[ind] + ownDt * (-0.33333333 * drho_dts[ind] * h[ind] / densities[ind]);
    float B_i = balsara[ind];
    float hi1 = 1 / h_i;
    float hi2x4 = h_i * h_i * 4;
    
    // Accumulators
    float3 dv_dt = 0;
    
    // Get our surrounding cells
    CellToScanRange range = setCellsToScanDynamic(positions[ind], ownDt * ownDt * length(a0) + h_i);
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
                if (distToCell2 > hi2x4) {
                  continue;
                }
                
                uint largeBoundCellIndex = cellPositionToLargeBoundIndex({ x, y, z });
                uint end = cellEnds[cellindex] + 1;
                // For each particle within the cell
                for (uint k = start; k < end; k++) {
                    const uint j = cellData[k].y;
                    if (largeBoundCellIndex != largeBoundCellParticles[j]) {
                        continue;
                    }

                    float3 a0 = accelerations[j] + grav_accelerations[j];
                    float3 x_j = positions[j] + ownDt * velocities[j] + ownDt * ownDt * a0;
                    
                    float3 x_ij = x_i - x_j;
                    float r = length(x_ij);
                    float h_j = h[j] + ownDt * (-0.33333333 * drho_dts[j] * h[j] / densities[j]);
                    if ((r > h_i + h_i and r > h_j + h_j) or r < 1e-12) {
                        // Cell is out of range or on top of us, skip.
                        continue;
                    }
                    // Fetch and calculate info about particle
                    float r1 = 1 / r;
                    float m_j = masses[j];
                    float rho_j = densities[j] + ownDt * drho_dts[j];
                    float u_j = internalEnergies[j] + ownDt * dInternalEnergy[j];
                    float4 pc = eos(materialIds[j], u_j, rho_j, texIron, texForesite, texFe, texHHe, texIce, texRock);
                    float fact_j = gradientTerms[j] * pc.x / (rho_j * rho_j);
                    float3 v_j = velocities[j] + ownDt * a0;
                    float c_j = pc.y;
                    float hj1 = 1 / h_j;
                    float3 v_ij = v_i - v_j;
                    float rho_ij = 0.5 * (rho_i + rho_j);
                    const float v_dot_r = dot(v_ij, x_ij);
                    float mu_ij = v_dot_r < 0 ? v_dot_r * r1 : 0;
                    float v_sigij = c_i + c_j - VISCOSITY_BETA * mu_ij;
                    
                    // Viscosity term
                    float B_j = balsara[j];
                    float B_ij = 0.5 * (B_i + B_j);
                    float nu_ij = - 0.5 * VISCOSITY_ALPHA * B_ij * mu_ij * v_sigij / rho_ij;
                    
                    float3 gW_i = gradW(x_ij, r, r1, hi1);
                    float3 gW_j = gradW(x_ij, r, r1, hj1);
                    float3 gW_ij = 0.5 * (gW_i + gW_j);
                    
                    // Acceleration terms
                    float3 t_i = fact_i * gW_i;
                    float3 t_j = fact_j * gW_j;
                    float3 t_visc = nu_ij * gW_ij;
                    
                    // Accumulate
                    dv_dt -= m_j * (t_i + t_j + t_visc);
                }
            }
        }
    }

    accelerations1[ind] = 0.5 * (accelerations[ind] + dv_dt);
}


