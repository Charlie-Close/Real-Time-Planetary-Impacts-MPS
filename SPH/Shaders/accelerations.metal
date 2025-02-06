//
//  accelerations.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include "utils/kernals/kernels.h"
#include "utils/cellsToScan/cellsToScan.h"
#include "utils/poles/poles.h"
#include "utils/morton.h"

constant float alpha = 1.5f;
constant float C_cfl = .7f;

//static inline float nu_ij(float3 x_ij, float3 v_ij, float h, float c, float rho) {
//    const float epsilon = 0.001;
//    const float v_dot_r = -1 * dot(v_ij, x_ij);
//    if (v_dot_r > 0) {
//        return 0;
//    }
//
//    const float mu = h * v_dot_r / (dot(x_ij, x_ij) + epsilon * h * h);
//    return (-alpha * c * mu + beta * mu * mu) / rho;
//}

constant float G = 6.67e-5;

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
                         device float* phiGrad,
                         device float* balsara,
                         device uint2* cellData,
                         device int* cellStarts,
                         device int* cellEnds,
                         device float* cellSize,
                         device int* cellsPerDim,
                         device bool* isAlive,
                         device atomic_int* dt,
                         uint index [[thread_position_in_grid]])
{
    int ind = cellData[index].y;
    if (!isAlive[ind]) {
        return;
    }
    float f_i = gradientTerms[ind];
    float P_i = pressures[ind];
    float rho_i = densities[ind];
    float fact_i = (f_i * P_i / (rho_i * rho_i));
    float3 x_i = positions[ind];
    float3 v_i = velocities[ind];
    float c_i = speedsOfSound[ind];
    float h_i = h[ind];
    float grav_f = f_i * phiGrad[ind];
    float B_i = balsara[ind];

    float3 dv_dt = 0;
    float du_dt = 0;
    float gradV = 0;
    
    int cellsPerDimT = *cellsPerDim;
    CellToScanRange range = setCellsToScanDynamic(x_i, *cellSize, *cellsPerDim, h_i);
    float v_sigi = 0.0001;
    
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
                    if (!isAlive[j]) {
                        continue;
                    }
                    float m_j = masses[j];
                    float f_j = gradientTerms[j];
                    float P_j = pressures[j];
                    float rho_j = densities[j];
                    float fact_j = (f_j * P_j / (rho_j * rho_j));
                    float3 x_j = positions[j];
                    float3 v_j = velocities[j];
                    float c_j = speedsOfSound[j];
                    float h_j = h[j];
//                    float grav_fj = f_j * phiGrad[j];
                    
                    float3 x_ij = x_i - x_j;
                    float l = length(x_ij);
                    if (l > h_i * 2 and l > h_j * 2) {
                        continue;
                    }
                    float3 v_ij = v_i - v_j;
                    float rho_ij = 0.5 * (rho_i + rho_j);
                    
                    const float v_dot_r = dot(v_ij, x_ij);
                    float mu_ij = v_dot_r < 0 ? v_dot_r / length(x_ij) : 0;
                    float v_sigij = c_i + c_j - 3 * mu_ij;
                    v_sigi = max(v_sigij, v_sigi);
                    
                    
                    float B_j = balsara[j];
                    float B_ij = 0.5 * (B_i + B_j);
                    float nu_ij = 0.5 * alpha * B_ij * mu_ij * v_sigij / rho_ij;
                    
                    float3 gW_i = - gradW(x_ij, h_i);
                    float3 gW_j = - gradW(x_ij, h_j);
                    
                    
                    float3 t_i = fact_i * gW_i;
                    float3 t_j = fact_j * gW_j;
                    float3 t_visc = nu_ij * 0.5 * (gW_i + gW_j);
                    
                    dv_dt -= m_j * (t_i + t_j + t_visc);
                    du_dt += m_j * (fact_i * dot(v_ij, gW_i) + 0.5 * nu_ij * dot(v_ij, 0.5 * (gW_i + gW_j)));
                }
            }
        }
    }

    accelerations[ind] += dv_dt;
    dInternalEnergy[ind] = du_dt;
    
    float goaldt = 2 * C_cfl * h_i / v_sigi;
    atomic_fetch_min_explicit(&(*dt), (int)ceil(10000 * goaldt), memory_order_relaxed);
}
