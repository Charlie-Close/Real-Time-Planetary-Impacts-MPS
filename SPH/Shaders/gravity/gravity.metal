//
//  gravity.metal
//  SPH
//
//  Created by Charlie Close on 23/01/2025.
//

#include "../utils/kernels.h"
#include "./poles/poles.h"
#include "helpers.h"

// Gravity is done in two stages, an up pass and a down pass on an octree. The octree is how we have
// grouped our particles. The root contains all particles. This then subdivides the particles into
// 8 branches. These branches subdivide again until a limit is hit and we have a leaf node, which just
// contains particles, and does not branch again.
//
// We first traverse up the octree from the lowest nodes (i.e. leaf nodes). Each node calculates it's
// multipole expansion. A first order multipole expansion is just the mass and center of mass. Higher
// orders allow more accurate approximations. The leaf nodes sum over their particles to produce an
// expansion (P2M - particle to multipole). The branch nodes first shift their child nodes multipoles
// to their center of mass, and then sums their expansion terms (M2M - multipole to multipole). The
// result of this is for every node in our tree we have a multipole which includes contriubtions from
// every node in the tree.
//
// We then traverse down the octree from the root. Each nodes checks if all other nodes are far enough
// away that we can use the multipole approximation. Then then produce a local expansion representing
// the resulting gravitational field at that point. A first order local field will just include the
// gravitational potential and the acceleration at the COM. Higher order expansion give more detail
// to the spacial derivative of this acceleration, giving greater accuracy. Each brach then shifts
// it's local expansion to each of its children's COM, and gives it to them to continue summing. We
// also store an array of all nodes which we were unable to use a multipole approximation on, so our
// children know which nodes to check (i.e. stops double counting). When we hit a leaf, we ensure there
// are no nodes whose gravitational fields have not been accounted for, and then we use our local
// expansion to caculate the acceleration due to gravity of each particle.

kernel void upPass(device float3* positions,
                   device float* masses,
                   device int* treeStructure,
                   device Multipole* multipoles,
                   device Local* locals,
                   device int* pointers,
                   device uint* parentIndexes,
                   device float* gravNorm,
                   device bool* active,
                   device int* nextActiveTime,
                   device int* globalTime,
                   uint index [[thread_position_in_grid]])
{
    int treePointer = pointers[index];
    int nParticles = treeStructure[treePointer];
    int dataPointer = treeStructure[treePointer + 1];
    
    if (nParticles == 0) {
        // We are looking at a branch node, so we sum the child multipoles (known in literature as
        // M2M - Multipole to Multipole)
        multipoles[dataPointer] = M2M(treeStructure, multipoles, active, parentIndexes, index, treePointer);
    } else {
        // We are looking at a leaf node, so we sum the child particles (known in literature as P2M
        // - Particle to Multipole)
        multipoles[dataPointer] = P2M(treeStructure, masses, positions, gravNorm, active, nextActiveTime, globalTime, treePointer);
    }
}

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
                     device bool* active,
                     uint index [[thread_position_in_grid]])
{
    int treePointer = pointers[index];
    
    // We are looking at the root. We won't even try to do a scan all the way up here.
    // We just zero out our children't local expansions.
    if (treePointer == 0) {
        zeroFirstBranches(treeStructure, multipoles, locals, unscannedIndexesOut);
        return;
    }
    
    
    int dataPointer = treeStructure[treePointer + 1];
    Multipole mp = multipoles[dataPointer];
    if (!mp.active) {
        // If we are inactive, no point calculating our local field.
        return;
    }
    
    // Get where we are going to store our unchecked nodes. Everytime we write to this array, we
    // increment our output pointer.
    uint outputPointer = index * MAX_UNCHECKED_POINTERS;
    uint maxOutputPointer = (index + 1) * MAX_UNCHECKED_POINTERS;
    
    int nParticles = treeStructure[treePointer];
    // We obviously wont check ourselves, however our children will need to check each other, so we need
    // to ourselves in the unchecked array.
    if (nParticles == 0) {
        unscannedIndexesOut[outputPointer] = treePointer;
        outputPointer++;
    }
    
    // Where our parent will have stored it's unchecked nodes
    uint parentIndex = parentIndexes[dataPointer];
    uint startScan = parentIndex * MAX_UNCHECKED_POINTERS;
    uint endScan = (parentIndex + 1) * MAX_UNCHECKED_POINTERS;
    
    // Stage 1: calculate our local expansion.
    Local local = locals[dataPointer];
    for (uint i = startScan; i < endScan; i++) {
        int toScan = unscannedIndexesIn[i];
        // If is -1, we have scanned everything we need to.
        if (toScan == -1) {
            break;
        }

        if (nParticles != 0) {
            // If we are looking at a leaf, we can't write to unscannedIndexes, so recursion is needed.
            recursiveScan(treeStructure, multipoles, mp, local, toScan, treePointer);
        } else {
            nonRecursiveScan(treeStructure, unscannedIndexesOut, multipoles, mp, local, outputPointer, maxOutputPointer, toScan, treePointer);
        }
    }
    
    if (outputPointer < maxOutputPointer and nParticles == 0) {
        // Set final index to -1 so children know where to stop scanning.
        unscannedIndexesOut[outputPointer] = -1;
    }
    
    
    
    // Stage 2: pass our local expansion down.
    if (nParticles == 0) {
        L2L(treeStructure, multipoles, locals, mp, local, treePointer);
    } else {
        localToLeaf(treeStructure, active, positions, accelerations, masses, gravNorm, multipoles, mp, local, treePointer, nParticles);
    }
}
