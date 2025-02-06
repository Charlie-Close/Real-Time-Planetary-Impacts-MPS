//
//  tree.metal
//  SPH
//
//  Created by Charlie Close on 23/01/2025.
//

#include "utils/kernals/kernels.h"
#include "utils/poles/poles.h"

constant float theta_lim = .7f;
constant float G = 6.67e-5;

kernel void upPass(device float3* positions,
                   device float* masses,
                   device int* treeStructure,
                   device Multipole* multipoles,
                   device Local* locals,
                   device int* pointers,
                   device uint* parentIndexes,
                   device float* gravNorm,
                   uint index [[thread_position_in_grid]])
{
    int treePointer = pointers[index];
    int nParticles = treeStructure[treePointer];
    int dataPointer = treeStructure[treePointer + 1];
    
    if (nParticles == 0) {
        Multipole mp;
        mp.minGrav = MAXFLOAT;
        mp.pos = { 0, 0, 0 };
        for (uint i = 0; i < N_EXPANSION_TERMS; i++) {
            mp.expansion[i] = 0;
        }
        float mass = 0;
        bool first = true;
        int start = treePointer + 2;
        int end = start + 8;
        for (int i = start; i < end; i++) {
            int childPointer = treeStructure[i];
            if (childPointer == -1) {
                continue;
            }
            int childDataPointer = treeStructure[childPointer + 1];
            parentIndexes[childDataPointer] = index;
            Multipole childMp = multipoles[childDataPointer];
            mp.minGrav = min(childMp.minGrav, mp.minGrav);
            mass += childMp.expansion[M];
            mp.pos += childMp.pos * childMp.expansion[M];
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
        
        float3 dims = mp.max - mp.min;
        mp.size = max3(dims.x, dims.y, dims.z);
        mp.pos /= mass;
        
        for (int i = start; i < end; i++) {
            int childPointer = treeStructure[i];
            if (childPointer == -1) {
                continue;
            }
            int childDataPointer = treeStructure[childPointer + 1];
            parentIndexes[childDataPointer] = index;
            Multipole childMp = multipoles[childDataPointer];
            float3 r = mp.pos - childMp.pos;
            Multipole transformed = transformMultipole(childMp, r);
            for (uint j = 0; j < N_EXPANSION_TERMS; j++) {
                mp.expansion[j] += transformed.expansion[j];
            }
        }
        addPowers(mp);
        
        multipoles[dataPointer] = mp;
    } else {
        multipoles[dataPointer] = P2M(treeStructure, masses, positions, gravNorm, treePointer);
    }
}


constant uint nUnscanned = 8192;

kernel void downPass(device float3* positions,
                     device float3* accelerations,
                     device float* masses,
                     device int* treeStructure,
                     device Multipole* multipoles,
                     device Local* locals,
                     device int* pointers,
                     device uint* parentIndexes,
                     device int* unscannedIndexesIn,
                     device int* unscannedIndexesOut,
                     device float* gravNorm,
                     device int* nextActiveTime,
                     device int* globalTime,
                     device bool* active,
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
        Local local;
        for (uint i = 0; i < N_EXPANSION_TERMS; i++) {
            local.expansion[i] = 0;
        }
        for (int i = 2; i < 10; i++) {
            int nodePointer = treeStructure[i];
            if (nodePointer == -1) {
                continue;
            }
            int nodeDataPointer = treeStructure[nodePointer + 1];
            float3 nodePos = multipoles[nodeDataPointer].pos;
            locals[nodeDataPointer].pos = nodePos;
            for (uint i = 0; i < N_EXPANSION_TERMS; i++) {
                locals[nodeDataPointer].expansion[i] = 0;
            }
        }
        unscannedIndexesOut[outputPointer] = -1;
        return;
    }
    
    
    
    int dataPointer = treeStructure[treePointer + 1];
    uint parentIndex = parentIndexes[dataPointer];
    uint startScan = parentIndex * nUnscanned;
    uint endScan = (parentIndex + 1) * nUnscanned;
    
    Multipole mp = multipoles[dataPointer];
    Local local = locals[dataPointer];
    
    for (uint i = startScan; i < endScan; i++) {
        int toScan = unscannedIndexesIn[i];
        // Check if we are at the end of the array.
        if (toScan == -1) {
            break;
        }
        
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
                Multipole nodeMp = multipoles[nodeDataPointer];
//                float theta = calculateAngularSize(mp.pos, multipoles, nodeDataPointer, mp.size);
                if (gravity_M2L_accept(mp, nodeMp) or nParticles != 0 or outputPointer == maxOutputPointer) {
//                if (false or nParticles != 0 or outputPointer == maxOutputPointer) {
                    Local nodeLocal = M2L(local.pos, nodeMp);
                    for (uint j = 0; j < N_EXPANSION_TERMS; j++) {
                        local.expansion[j] += nodeLocal.expansion[j];
                    }
                } else {
                    unscannedIndexesOut[outputPointer] = nodePointer;
                    outputPointer++;
                }
            }
        } else {
            // We are looking at a leaf.
            int nodeDataPointer = treeStructure[toScan + 1];
            Multipole nodeMp = multipoles[nodeDataPointer];
            Local nodeLocal = M2L(local.pos, nodeMp);
            for (uint j = 0; j < N_EXPANSION_TERMS; j++) {
                local.expansion[j] += nodeLocal.expansion[j];
            }
        }
    }
    
    
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
            float3 r = multipoles[nodeDataPointer].pos - local.pos;
            locals[nodeDataPointer] = transformLocal(local, r);
        }
    } else {
        // Give acceleration to all child particles:
        int start = treePointer + 2;
        int end = start + nParticles;
        for (int i = start; i < end; i++) {
            int particlePointer = treeStructure[i];
            if (!active[particlePointer]) {
                continue;
            }
            float3 x_i = positions[particlePointer];
            float3 particleAcceleration = L2P(local, x_i);
                        
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
                particleAcceleration -= (m_j * dphi_dr(x_ij, h_i)) * r_hat;
            }
            accelerations[particlePointer] = G * particleAcceleration;
            gravNorm[particlePointer] = length(particleAcceleration);
        }
    }
}
