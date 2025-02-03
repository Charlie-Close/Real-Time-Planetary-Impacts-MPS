//
//  tree.metal
//  SPH
//
//  Created by Charlie Close on 23/01/2025.
//

#include "utils/kernals/kernels.h"
#include "utils/poles/poles.h"

constant float theta_lim = .3f;
constant float G = 6.67e-5;

static inline float calculateAngularSize(float3 x_i, device Multipole* treeData, int pointer, int ownSize) {
    float size = treeData[pointer].size;
    float3 x_j = treeData[pointer].pos;
    return size / max(fast::length(x_i - x_j) - 0.5 * ownSize, 1e-12f);
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
    const float3 r_hat = normalize(x_ij);
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

        const float3 r_hat = normalize(x_ij);
        acceleration -= (G * m_j * dphi_dr(x_ij, h_i)) * r_hat;
        phiGrad += m_j * dphi_dh(x_ij, h_i);
    }
}

kernel void upPass(device float3* positions,
                   device float* masses,
                   device int* treeStructure,
                   device Multipole* treeData,
                   device int* pointers,
                   device uint* parentIndexes,
                   uint index [[thread_position_in_grid]])
{
    int treePointer = pointers[index];
    int nParticles = treeStructure[treePointer];
    int dataPointer = treeStructure[treePointer + 1];
    Multipole mp;
    mp.mass = 0;
    mp.pos = 0;
    bool first = true;
    
    if (nParticles == 0) {
        int start = treePointer + 2;
        int end = start + 8;
        for (int i = start; i < end; i++) {
            int childPointer = treeStructure[i];
            if (childPointer == -1) {
                continue;
            }
            int childDataPointer = treeStructure[childPointer + 1];
            parentIndexes[childDataPointer] = index;
            Multipole childMp = treeData[childDataPointer];
            mp.mass += childMp.mass;
            mp.pos += childMp.pos * childMp.mass;
            if (first) {
                mp.max = childMp.max;
                mp.min = childMp.min;
                first = false;
            } else {
                mp.max.x = max(childMp.max.x, mp.max.x);
                mp.max.y = max(childMp.max.y, mp.max.y);
                mp.max.z = max(childMp.max.z, mp.max.z);
                mp.min.x = min(childMp.min.x, mp.min.x);
                mp.min.y = min(childMp.min.y, mp.min.y);
                mp.min.z = min(childMp.min.z, mp.min.z);
            }
        }
    } else {
        // We are looking at a leaf
        int start = treePointer + 2;
        int end = start + nParticles;
        for (int i = start; i < end; i++) {
            int p = treeStructure[i];
            float p_mass = masses[p];
            float3 p_position = positions[p];
            if (first) {
                mp.max = p_position;
                mp.min = p_position;
                first = false;
            } else {
                mp.max.x = max(p_position.x, mp.max.x);
                mp.max.y = max(p_position.y, mp.max.y);
                mp.max.z = max(p_position.z, mp.max.z);
                mp.min.x = min(p_position.x, mp.min.x);
                mp.min.y = min(p_position.y, mp.min.y);
                mp.min.z = min(p_position.z, mp.min.z);
            }
            
            
            mp.mass += p_mass;
            mp.pos += p_position * p_mass;
        }
    }
    
    mp.pos /= mp.mass;
    float3 dims = mp.max - mp.min;
    mp.size = max3(dims.x, dims.y, dims.z);

    treeData[dataPointer] = mp;
}







constant uint nUnscanned = 8192;

kernel void downPass(device float3* positions,
                     device float3* accelerations,
                     device float* masses,
                     device int* treeStructure,
                     device Multipole* treeData,
                     device int* pointers,
                     device uint* parentIndexes,
                     device int* unscannedIndexesIn,
                     device int* unscannedIndexesOut,
                     device atomic_int* maxCount,
                     device bool* finalLayer,
                     uint index [[thread_position_in_grid]])
{
    int treePointer = pointers[index];
    uint outputPointer = index * nUnscanned;
    uint maxOutputPointer = (index + 1) * nUnscanned;
    int nParticles = treeStructure[treePointer];
    
    // We won't scan ourselves, so get children to do it.
    if (nParticles == 0) {
        unscannedIndexesOut[outputPointer] = treePointer;
        outputPointer++;
    }
    
    if (treePointer == 0) {
        // We are looking at the root. We won't even try to do a scan all the way up here.
        for (int i = 2; i < 10; i++) {
            int nodePointer = treeStructure[i];
            if (nodePointer == -1) {
                continue;
            }
            int nodeDataPointer = treeStructure[nodePointer + 1];
            treeData[nodeDataPointer].external = { 0, 0, 0 };
        }
        unscannedIndexesOut[outputPointer] = -1;
        return;
    }
    
    
    
    int dataPointer = treeStructure[treePointer + 1];
    uint parentIndex = parentIndexes[dataPointer];
    uint startScan = parentIndex * nUnscanned;
    uint endScan = (parentIndex + 1) * nUnscanned;
    
    Multipole mp = treeData[dataPointer];
    float3 acceleration = treeData[dataPointer].external;
    float phiGrad = 0;
    int nToScan = 0;
    
    
    
    
    for (uint i = startScan; i < endScan; i++) {
        int toScan = unscannedIndexesIn[i];
        // Check if we are at the end of the array.
        if (toScan == -1) {
            break;
        }
        nToScan++;
        
        if (treeStructure[toScan] == 0) {
            // We are looking at a branch.
            // Check the 8 children.
            int start = toScan + 2;
            int end = start + 8;
            for (int j = start; j < end; j++) {
                int nodePointer = treeStructure[j];
                if (nodePointer == treePointer or nodePointer == -1) {
                    // Ignore ourselves and empty nodes
                    continue;
                }
                int nodeDataPointer = treeStructure[nodePointer + 1];
                float theta = calculateAngularSize(mp.pos, treeData, nodeDataPointer, mp.size);
                if (theta < theta_lim or nParticles != 0 or outputPointer == maxOutputPointer) {
                    getAccelerationFromBranch(acceleration, phiGrad, 0.05, treeData, nodeDataPointer, mp.pos);
                } else {
                    unscannedIndexesOut[outputPointer] = nodePointer;
                    outputPointer++;
                }
            }
        } else {
            // We are looking at a leaf.
            // Loop through all the particles in the leaf.
            sumAccelerationsInLeaf(acceleration, phiGrad, positions, masses, treeStructure, toScan, mp.pos, 0.05);
        }
    }
    
//    atomic_fetch_max_explicit(&(*maxCount), nToScan, memory_order_relaxed);
    
    if (outputPointer < maxOutputPointer and nParticles == 0) {
        unscannedIndexesOut[outputPointer] = -1;
    }
    
    if (nParticles == 0) {
        // Give acceleration to all child branches
        int start = treePointer + 2;
        int end = start + 8;
        for (int j = start; j < end; j++) {
            int nodePointer = treeStructure[j];
            int nodeDataPointer = treeStructure[nodePointer + 1];
            treeData[nodeDataPointer].external = acceleration;
        }
    } else {
        // Give acceleration to all child particles:
        int start = treePointer + 2;
        int end = start + nParticles;
        for (int i = start; i < end; i++) {
            int particlePointer = treeStructure[i];
            float3 particleAcceleration = { acceleration.x, acceleration.y, acceleration.z };
            float3 x_i = positions[particlePointer];
            float h_i = 0.05;
            for (int j = start; j < end; j++) {
                if (i == j) {
                    continue;
                }
                int p_j = treeStructure[j];
                float3 x_j = positions[p_j];
                float m_j = masses[p_j];

                float3 x_ij = x_i - x_j;
                const float r = fast::length(x_ij);
                if (r == 0) {
                    continue;
                }

                const float3 r_hat = x_ij / r;
                particleAcceleration -= (G * m_j * dphi_dr(x_ij, h_i)) * r_hat;
            }
            accelerations[particlePointer] = particleAcceleration;
        }
    }
}
