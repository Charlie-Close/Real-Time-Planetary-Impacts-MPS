//
//  accelerations.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include "utils/kernels.h"
#include "utils/cellsToScan.h"
#include "utils/morton.h"
#include "../Parameters.h"

kernel void acceleration(device float3* positions,
                         device float3* velocities,
                         device float3* accelerations,
                         device float* densities,
                         device float* internalEnergies,
                         device float* masses,
                         device float* pressures,
                         device float* h,
                         device float* gradientTerms,
                         device float* speedsOfSound,
                         device float* dInternalEnergy,
                         device float* balsara,
                         device float* dhdts,
                         device uint2* cellData,
                         device int* cellStarts,
                         device int* cellEnds,
                         device float* cellSize,
                         device int* cellsPerDim,
                         device bool* active,
                         device int* nextActiveTime,
                         device int* globalTime,
                         device atomic_int* dt,
                         uint index [[thread_position_in_grid]])
{
    int ind = cellData[index].y;
    // If we are inactive this timestep, we skip
    if (!active[ind]) {
        return;
    }
    // Get all of our relevant terms
    float f_i = gradientTerms[ind];
    float P_i = pressures[ind];
    float rho_i = densities[ind];
    float fact_i = (f_i * P_i / (rho_i * rho_i));
    float3 x_i = positions[ind];
    float3 v_i = velocities[ind];
    float c_i = speedsOfSound[ind];
    float h_i = h[ind];
    float B_i = balsara[ind];
    float hi1 = 1 / h_i;
    float hi2x4 = h_i * h_i * 4;

    // Accumulators
    float3 dv_dt = 0;
    float du_dt = 0;
    float dhdt = 0;
    
    int cellsPerDimT = *cellsPerDim;
    // Get our surrounding cells
    CellToScanRange range = setCellsToScanDynamic(x_i, *cellSize, *cellsPerDim, h_i);
    float v_sigi = 1e-12;
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
                uint end = cellEnds[cellindex] + 1;
                // For each particle within the cell
                for (uint k = start; k < end; k++) {
                    const uint j = cellData[k].y;
                    float3 x_j = positions[j];
                    float3 x_ij = x_i - x_j;
                    float r2 = length_squared(x_ij);
                    float h_j = h[j];
                    float hj2x4 = h_j * h_j * 4;
                    if (r2 > hi2x4 and r2 > hj2x4) {
                        // Cell is out of range, skip.
                        continue;
                    }
                    // Fetch and calculate info about particle
                    float r = sqrt(r2);
                    float r1 = 1 / r;
                    float m_j = masses[j];
                    float f_j = gradientTerms[j];
                    float P_j = pressures[j];
                    float rho_j = densities[j];
                    float fact_j = (f_j * P_j / (rho_j * rho_j));
                    float3 v_j = velocities[j];
                    float c_j = speedsOfSound[j];
                    float hj1 = 1 / h_j;
                    float3 v_ij = v_i - v_j;
                    float rho_ij = 0.5 * (rho_i + rho_j);
                    const float v_dot_r = dot(v_ij, x_ij);
                    float mu_ij = v_dot_r < 0 ? v_dot_r * r1 : 0;
                    float v_sigij = c_i + c_j - 3 * mu_ij;
                    v_sigi = max(v_sigij, v_sigi);
                    
                    // Viscosity term
                    float B_j = balsara[j];
                    float B_ij = 0.5 * (B_i + B_j);
                    float nu_ij = 0.5 * VISCOSITY_ALPHA * B_ij * mu_ij * v_sigij / rho_ij;
                    
                    float3 gW_i = - gradW(x_ij, r, r1, hi1);
                    float3 gW_j = - gradW(x_ij, r, r1, hj1);
                    
                    // Acceleration terms
                    float3 t_i = fact_i * gW_i;
                    float3 t_j = fact_j * gW_j;
                    float3 t_visc = nu_ij * 0.5 * (gW_i + gW_j);
                    
                    // Accumulate
                    dv_dt -= m_j * (t_i + t_j + t_visc);
                    du_dt += m_j * (fact_i * dot(v_ij, gW_i) + 0.5 * nu_ij * dot(v_ij, 0.5 * (gW_i + gW_j)));
                    dhdt += (m_j / rho_j) * dot(v_ij, gW_i);
                }
            }
        }
    }

    accelerations[ind] += dv_dt;
    dInternalEnergy[ind] = du_dt;
    dhdts[ind] = dhdt;

    // Time stepping criterion.
    float goaldt = 2.0f * CFL * h_i / v_sigi;
    int r = 0;
    float dtCandidate = DT0;
    while (dtCandidate > goaldt && r < R_MAX) {
        dtCandidate *= 0.5f;
        r++;
    }
    int integerDt = (int)ceil(DT_MIN * dtCandidate);
    nextActiveTime[ind] = (*globalTime) + integerDt;
    // dt is an integer as atomic min doesn't seem to like floats.
    atomic_fetch_min_explicit(&(*dt), integerDt, memory_order_relaxed);
}
