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
                         device float* dh_dts,
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
                         device float* alpha,
                         device float* da_dt,
                         device float* alphaLoc,
                         device float* pAlphaLoc,
                         uint index [[thread_position_in_grid]],
                         uint localId [[thread_position_in_threadgroup]],
                         uint threadsPerGroup [[threads_per_threadgroup]])
{
    uint ind = cellData[index].y;
    float3 x_i = positions[ind];
    
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
    float hSupportSquared = GAMMA * GAMMA * h_i * h_i;
    float alpha_i = alpha[ind];
    float H1 = hi1 * gamma1;
    float W_kernal_fac = KERNEL_CONSTANT * H1 * H1 * H1;
    float gW_kernal_fac = KERNEL_CONSTANT * H1 * H1 * H1 * H1;
    
    

    // Accumulators
    float3 dv_dt(0);
    float du_dt = 0;
    float dh_dt = 0;
    float alpha_loc = (masses[ind] / rho_i) * pAlphaLoc[ind] * W_fast(0, H1, W_kernal_fac);
    
    // Get our surrounding cells
    CellToScanRange range = setCellsToScanDynamic(x_i, min(ACC_LOOK_EXTEND * h_i, MAX_SMOOTHING_LENGTH));
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
                if (distToCell2 > ACC_LOOK_EXTEND * ACC_LOOK_EXTEND * hSupportSquared) {
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

                    float3 x_j = positions[j];
                    
                    float3 x_ij = x_i - x_j;
                    float r2 = length_squared(x_ij);
                    float h_j = h[j];
                    float hjSupport = GAMMA * h_j;
                    if ((r2 > hSupportSquared and r2 > hjSupport * hjSupport) or r2 < 1e-12) {
                        // Cell is out of range or on top of us, skip.
                        continue;
                    }
                    // Fetch and calculate info about particle
                    float r = sqrt(r2);
                    float r1 = 1 / r;
                    float m_j = masses[j];
                    float rho_j = densities[j];
                    float m_rho = m_j / rho_j;
                    float fact_j = pressures[j] * gradientTerms[j] / (rho_j * rho_j);
                    float3 v_ij = v_i - velocities[j];
                    const float v_dot_r = dot(v_ij, x_ij);
                    float mu_ij = v_dot_r < 0 ? v_dot_r * r1 : 0;
                    float v_sigij = c_i + speedsOfSound[j] - VISCOSITY_BETA * mu_ij;
                    v_sigi = max(v_sigij, v_sigi);
                    float nu_ij = - 0.25 * (alpha_i + alpha[j]) * (B_i + balsara[j]) * mu_ij * v_sigij / (rho_i + rho_j);

                    
                    float3 gW_i = gradW_fast(x_ij, r, r1, H1, gW_kernal_fac);
                    float3 gW_j = gradW(x_ij, r, r1, 1 / h_j);
                    float W_i = W_fast(r, H1, W_kernal_fac);
                    
                    // Acceleration terms
                    const float3 visc_acc_term = 0.5f * nu_ij * (gW_i + gW_j);
                    const float3 SPH_acc_term = (fact_i * gW_i + fact_j * gW_j);
                    
                    // Accumulate
                    dv_dt -= m_j * (SPH_acc_term + visc_acc_term);
                    du_dt += m_j * (fact_i * dot(v_ij, gW_i) + 0.5 * dot(v_ij, visc_acc_term));
                    dh_dt += m_rho * dot(v_ij, gW_i);
                    alpha_loc += m_rho * pAlphaLoc[j] * W_i;
                }
            }
        }
    }
    
    const float tau = h_i / (0.1 * v_sigi);
    alpha_loc = clamp(alpha_loc, ALPHA_MIN, ALPHA_MAX);
    if (alpha_i < alpha_loc) {
        alpha_i = alpha_loc;
    }
    alpha[ind] = alpha_i;
    da_dt[ind] = (alpha_loc - alpha_i) / tau;
    alphaLoc[ind] = alpha_loc;

    // Time stepping criterion.
    float dtCFL = 2.0f * CFL * GAMMA * h_i / v_sigi;
    float goaldt = clamp(dtCFL, MIN_DT, MAX_DT);
    int integerDt = 1;
    while (2 * MIN_DT * integerDt < goaldt and 2 * MIN_DT * integerDt < MAX_DT) {
        integerDt += integerDt;
    }
    
    nextActiveTime[ind] = globalTime + integerDt;
    // dt is an integer as atomic min doesn't work with floats.
    atomic_fetch_min_explicit(&(*dt), integerDt, memory_order_relaxed);
    dInternalEnergy[ind] = du_dt;
    dh_dts[ind] = - 0.333333333 * h_i * dh_dt;
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
                             device float* dh_dts,
                             device uint2* cellData,
                             device int* cellStarts,
                             device int* cellEnds,
                             device uint* largeBoundCellParticles,
                             device bool* active,
                             device bool* alive,
                             device int* nextActiveTime,
                             device int& globalTime,
                             device int& dt,
                             device float* pressures,
                             device int* materialIds,
                             device float* alpha,
                             device float* da_dt,
                             device float* alphaLoc,
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
    
    // Cap time step at 8x the smallest timestep
    float globalTimeStep = MIN_DT * dt;
    float ownDt = MIN_DT * (nextActiveTime[ind] - globalTime);
    ownDt = min(ownDt, 4 * globalTimeStep);
    nextActiveTime[ind] = globalTime + (int)(ownDt / MIN_DT);
    float half_dt = 0.5 * ownDt;
        
    // Get all of our relevant terms
    float3 a_i = accelerations[ind] + grav_accelerations[ind];
    float3 x_i = positions[ind] + half_dt * velocities[ind] + 0.5 * half_dt * half_dt * a_i;
    float rho_i = densities[ind] * exp(-(3 / h[ind]) * dh_dts[ind] * half_dt);
    float3 v_i = velocities[ind] + half_dt * a_i;
    float u = max(internalEnergies[ind] + half_dt * dInternalEnergy[ind], 0.f);
    float4 pc = eos(materialIds[ind], u, rho_i, texIron, texForesite, texFe, texHHe, texIce, texRock);
    float fact_i = gradientTerms[ind] * pc.x / (rho_i * rho_i);
    float c_i = pc.y;
    float h_i = clamp(h[ind] * exp((1 / h[ind]) * dh_dts[ind] * half_dt), 1e-12, MAX_SMOOTHING_LENGTH);
    float B_i = balsara[ind];
    float hi1 = 1 / h_i;
    float hSupportSquared = GAMMA * GAMMA * h_i * h_i;
    float alpha_i = clamp(alphaLoc[ind] + (alpha[ind] - alphaLoc[ind]) * exp(da_dt[ind]), ALPHA_MIN, ALPHA_MAX);
    float H1 = hi1 * gamma1;
    float gW_kernal_fac = KERNEL_CONSTANT * H1 * H1 * H1 * H1;
    
    // Accumulators
    float3 dv_dt = 0;
    float du_dt = 0;
    float dh_dt = 0;

    // Get our surrounding cells
    CellToScanRange range = setCellsToScanDynamic(positions[ind], half_dt * half_dt * length(a_i) + min(ACC_LOOK_EXTEND * h_i, MAX_SMOOTHING_LENGTH));
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
                if (distToCell2 > ACC_LOOK_EXTEND * ACC_LOOK_EXTEND * hSupportSquared) {
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

                    float3 a_j = accelerations[j] + grav_accelerations[j];
                    float3 x_j = positions[j] + half_dt * velocities[j] + 0.5 * half_dt * half_dt * a_j;
                    
                    float3 x_ij = x_i - x_j;
                    float r2 = length_squared(x_ij);
                    float h_j = h[j];
                    h_j = clamp(h_j * exp(dh_dts[j] * half_dt / h_j), 1e-12, MAX_SMOOTHING_LENGTH);
                    float hjSupport = GAMMA * h_j;
                    if ((r2 > hSupportSquared and r2 > hjSupport * hjSupport) or r2 < 1e-12) {
                        // Cell is out of range or on top of us, skip.
                        continue;
                    }
                    
                    // Fetch and calculate info about particle
                    float r = sqrt(r2);
                    float r1 = 1 / r;
                    float m_j = masses[j];
                    float rho_j = densities[j] * exp(-(3 / h[j]) * dh_dts[j] * half_dt);
                    float u_j = max(internalEnergies[j] + half_dt * dInternalEnergy[j], 0.f);
                    float4 pc = eos(materialIds[j], u_j, rho_j, texIron, texForesite, texFe, texHHe, texIce, texRock);
                    float fact_j = gradientTerms[j] * pc.x / (rho_j * rho_j);
                    float3 v_j = velocities[j] + half_dt * a_j;
                    float3 v_ij = v_i - v_j;
                    const float v_dot_r = dot(v_ij, x_ij);
                    float mu_ij = v_dot_r < 0 ? v_dot_r * r1 : 0;
                    float alpha_j = clamp(alphaLoc[j] + (alpha[j] - alphaLoc[j]) * exp(da_dt[j]), ALPHA_MIN, ALPHA_MAX);
                    float v_sigij = c_i + pc.y - VISCOSITY_BETA * mu_ij;
                    v_sigi = max(v_sigij, v_sigi);
                                    
                    float nu_ij = - 0.25 * (alpha_i + alpha_j) * (B_i + balsara[j]) * mu_ij * v_sigij / (rho_i + rho_j);
                    
                    float3 gW_i = gradW_fast(x_ij, r, r1, H1, gW_kernal_fac);
                    float3 gW_j = gradW(x_ij, r, r1, 1 / h_j);
                    
                    // Acceleration terms
                    const float3 visc_acc_term = 0.5f * nu_ij * (gW_i + gW_j);
                    const float3 SPH_acc_term = (fact_i * gW_i + fact_j * gW_j);
                    
                    // Accumulate
                    dv_dt -= m_j * (SPH_acc_term + visc_acc_term);
                    du_dt += m_j * (fact_i * dot(v_ij, gW_i) + 0.5 * dot(v_ij, visc_acc_term));
                    dh_dt += (m_j / rho_j) * dot(v_ij, gW_i);
                }
            }
        }
    }

    accelerations1[ind] = dv_dt;
    dInternalEnergy[ind] = du_dt;
    dh_dt = - 0.333333333 * h_i * dh_dt;
    // Step dh_dt back by dt
    float prev_h_pred = h_i * exp(- (1 / h_i) * dh_dt * half_dt);
    float prev_dh_dt = dh_dt * prev_h_pred / h_i;
    dh_dts[ind] = 0.5 * (dh_dts[ind] + prev_dh_dt);
}


