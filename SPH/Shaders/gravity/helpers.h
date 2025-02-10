//
//  helpers.h
//  SPH
//
//  Created by Charlie Close on 07/02/2025.
//

#ifndef helpers_h
#define helpers_h

#include <metal_stdlib>
#include "./poles/poles.h"
using namespace metal;

static inline void zeroFirstBranches(device int* treeStructure, device Multipole* multipoles, device Local* locals, device int* unscannedIndexesOut) {
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
    // Make children scan root, and then stop.
    unscannedIndexesOut[0] = 0;
    unscannedIndexesOut[1] = -1;
    return;
}

static inline bool checkAndAddLocalExpansion(device int* treeStructure, device Multipole* multipoles, thread Multipole& mp, thread Local& local, int toScan, bool forceAccept) {
    int nodeDataPointer = treeStructure[toScan + 1];
    Multipole nodeMp = multipoles[nodeDataPointer];
    if (forceAccept or gravity_M2L_accept(mp, nodeMp)) {
        Local nodeLocal = M2L(local.pos, nodeMp);
        for (uint j = 0; j < N_EXPANSION_TERMS; j++) {
            local.expansion[j] += nodeLocal.expansion[j];
        }
        return true;
    }
    return false;
}

static inline void localToLeaf(device int* treeStructure, device bool* active, device float3* positions, device float3* accelerations, device float* masses, device float* gravNorm, device Multipole* multipoles, thread Multipole& mp, thread Local& local, int treePointer, int nParticles) {
    // Give acceleration to all child particles:
    int start = treePointer + 2;
    int end = start + nParticles;
    for (int i = start; i < end; i++) {
        int particlePointer = treeStructure[i];
        if (!active[particlePointer]) {
            // If particle is inactive, we skip it.
            continue;
        }
        // Otherwise, sum over contributions from surrounding particles.
        float3 x_i = positions[particlePointer];
        float3 particleAcceleration = L2P(local, x_i);
                    
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
            particleAcceleration -= (m_j * dphi_dr(x_ij, GRAVITY_SMOOTHING_LENGTH)) * r_hat;
        }
        accelerations[particlePointer] = G * particleAcceleration;
        gravNorm[particlePointer] = length(particleAcceleration);
    }
}

static inline void recursiveScan(device int* treeStructure, device Multipole* multipoles, thread Multipole& mp, thread Local& local, int toScan, int treePointer) {
    // If we can treat the parent as a multipole, just do that and continue
    if (checkAndAddLocalExpansion(treeStructure, multipoles, mp, local, toScan, treeStructure[toScan] != 0)) {
        return;
    }
    
    int indexes[GRAVITY_MAX_RECURSION];
    int ends[GRAVITY_MAX_RECURSION];
    uint currentDepth = 0;
    indexes[0] = toScan + 2;
    ends[0] = indexes[0] + 8;
        
    int count = 0;
    while (indexes[0] < ends[0] and count < 2000) {
        count++;
        int nodePointer = treeStructure[indexes[currentDepth]];
        if (nodePointer == treePointer or nodePointer == -1) {
            // Ignore ourselves and empty nodes.
            indexes[currentDepth]++;
            while (indexes[currentDepth] == ends[currentDepth] and currentDepth != 0) {
                currentDepth--;
                indexes[currentDepth]++;
            }
            continue;
        }
        
        if (!checkAndAddLocalExpansion(treeStructure, multipoles, mp, local, nodePointer, treeStructure[nodePointer] != 0 or currentDepth == GRAVITY_MAX_RECURSION - 1)) {
//          If we can't treat as multipole, step down the tree.
            currentDepth++;
            indexes[currentDepth] = nodePointer + 2;
            ends[currentDepth] = indexes[currentDepth] + 8;
            continue;
        }
        
        indexes[currentDepth]++;
        while (indexes[currentDepth] == ends[currentDepth] and currentDepth != 0) {
            currentDepth--;
            indexes[currentDepth]++;
        }
    }
}

void nonRecursiveScan(device int* treeStructure, device int* unscannedIndexesOut, device Multipole* multipoles, thread Multipole& mp, thread Local& local, thread uint& outputPointer, uint maxOutputPointer, int toScan, int treePointer) {
    // If we can treat the parent as a multipole, just do that and continue
    if (checkAndAddLocalExpansion(treeStructure, multipoles, mp, local, toScan, false)) {
        return;
    }
    
    // We couldn't treat the entire node as a multipole, so we will check it's children.
    int start = toScan + 2;
    int end = start + 8;
    for (int j = start; j < end; j++) {
        int childPointer = treeStructure[j];
        if (childPointer == treePointer or childPointer == -1) {
            // Ignore ourselves and empty nodes
            continue;
        }
        
        if (!checkAndAddLocalExpansion(treeStructure, multipoles, mp, local, childPointer, treeStructure[childPointer] != 0 or outputPointer == maxOutputPointer)) {
            // If we can't treat as multipole, add to unchecked indexes.
            unscannedIndexesOut[outputPointer] = childPointer;
            outputPointer++;
        }
    }
}


#endif /* helpers_h */
