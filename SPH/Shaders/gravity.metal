//
//  gravity.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//
#include "utils/kernals/kernels.h"
#include "utils/poles/poles.h"

constant float theta_lim = .7f;
constant float G = 6.67e-5;

static inline float calculateAngularSize(float3 x_i, device Multipole* treeData, int pointer) {
    float size = treeData[pointer].size;
    float3 x_j = treeData[pointer].pos;
    return size / fast::length(x_i - x_j);
}

static inline void getAccelerationFromBranch(thread float3& acceleration,
                               thread float& phiGrad,
                               float h_i,
                               device Multipole* treeData,
                               int pointer,
                               float3 x_i)
{
    Multipole mp = treeData[pointer];
    
    float3 x_ij = x_i - mp.pos;
    const float r = fast::length(x_ij);
    
    if (r == 0) {
        return;
    }
    const float3 r_hat = x_ij / r;
    acceleration -= (G * mp.mass * dphi_dr(x_ij, h_i)) * r_hat;
    phiGrad += mp.mass * dphi_dh(x_ij, h_i);
}

static inline void sumAccelerationsInLeaf(thread float3& acceleration,
                            thread float& phiGrad,
                              device float3* positions,
                              device float* masses,
                              device int* treeStructure,
                              int leafPointer,
                              float3 x_i,
                            float h_i)
{
    int nParticles = treeStructure[leafPointer];
    int start = leafPointer + 2;
    int end = start + nParticles;
    
    for (int j = start; j < end; j++) {
        int p_j = treeStructure[j];
        float3 x_j = positions[p_j];
        float m_j = masses[p_j];

        float3 x_ij = x_i - x_j;
        const float r = fast::length(x_ij);
        if (r == 0) {
            continue;
        }

        const float3 r_hat = x_ij / r;
        acceleration -= (G * m_j * dphi_dr(x_ij, h_i)) * r_hat;
        phiGrad += m_j * dphi_dh(x_ij, h_i);
    }
}

constant int maxStackDepth = 12;

float4 findAccelerationFromRoot(device float3* positions,
                                device float* masses,
                                float h_i,
                                device int* treeStructure,
                                device Multipole* treeData,
                                float3 x_i)
{
    float3 acceleration = { 0, 0, 0 };
    float phiGrad = 0;
    int currentDepth = 0;
    
    int curIndex[maxStackDepth];
    int end[maxStackDepth];
    
    for (int i = 0; i < maxStackDepth; i++) {
        end[i] = 0;
        curIndex[i] = 0;
    }
    
    curIndex[0] = 2;
    end[0] = curIndex[0] + 8;
    
    while (curIndex[0] < end[0]) {
        // path[currentDepth] is node we are scanning. curIndex[currentDepth] is the index of the node we are at.
        int treePointer = treeStructure[curIndex[currentDepth]];
        if (treePointer == -1) {
            // Nothing here.
            curIndex[currentDepth] += 1;
            while (curIndex[currentDepth] >= end[currentDepth] and currentDepth > 0) {
                // If we are at the end, move back up the tree.
                currentDepth -= 1;
                curIndex[currentDepth] += 1;
            }
            continue;
        }

        if (treeStructure[treePointer] == 0) {
            // We are looking at a branch
            int dataPointer = treeStructure[treePointer + 1];
            float theta = calculateAngularSize(x_i, treeData, dataPointer);
            if (theta < theta_lim or currentDepth == maxStackDepth - 1) {
                getAccelerationFromBranch(acceleration, phiGrad, h_i, treeData, dataPointer, x_i);
                curIndex[currentDepth] += 1;
                while (curIndex[currentDepth] >= end[currentDepth] and currentDepth > 0) {
                    // If we are at the end, move back up the tree.
                    currentDepth -= 1;
                    curIndex[currentDepth] += 1;
                }
                continue;
            }
            // Step down the tree
            currentDepth += 1;
            curIndex[currentDepth] = treePointer + 2;
            end[currentDepth] = curIndex[currentDepth] + 8;
            continue;
        }

        // We are looking at a leaf.
        int dataPointer = treeStructure[treePointer + 1];
        float theta = calculateAngularSize(x_i, treeData, dataPointer);
        if (theta < theta_lim) {
            getAccelerationFromBranch(acceleration, phiGrad, h_i, treeData, dataPointer, x_i);
        } else {
            sumAccelerationsInLeaf(acceleration, phiGrad, positions, masses, treeStructure, treePointer, x_i, h_i);
        }
        
        curIndex[currentDepth] += 1;
        while (curIndex[currentDepth] >= end[currentDepth] and currentDepth > 0) {
            // If we are at the end, move back up the tree.
            currentDepth -= 1;
            curIndex[currentDepth] += 1;
        }
    }

    return float4(acceleration, phiGrad);
}



kernel void gravity(device float3* positions,
                     device float* masses,
                     device float3* accelerations,
                     device float* h,
                     device float* phiGrad,
                     device int* treeStructure,
                     device Multipole* poles,
                     device uint2* particleCells,
                     device bool* isAlive,
                     uint index [[thread_position_in_grid]])
{
    int ind = particleCells[index].y;
    if (!isAlive[ind]) {
        return;
    }
    float h_i = h[ind];
    float3 x_i = positions[ind];
    float4 gravityData = findAccelerationFromRoot(positions, masses, h_i, treeStructure, poles, x_i);
    accelerations[ind] = gravityData.xyz;
    phiGrad[ind] = gravityData.w;
}

