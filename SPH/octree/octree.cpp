//
//  octree.cpp
//  SPH
//
//  Created by Charlie Close on 03/11/2024.
//

#include "octree.hpp"
#include <limits>
#include <algorithm>
#include <iostream>

// A simple bounding box struct
struct BoundingBox {
    simd_float3 min;
    simd_float3 max;
};

/**
 * Recursively builds the octree for a subset of particles.
 *
 * @param positions       Pointer to all particle positions.
 * @param subsetIndices   Which particles (indices) belong in this subtree.
 * @param level           Current depth in the tree.
 * @param maxLeafSize     Max #particles allowed in a leaf before splitting.
 * @param octreeData      Flat integer array describing the octree structure.
 * @param treeLevels      2D array: treeLevels[d] lists node indices at depth d.
 * @param nodeValues      Float array for branch node data (7 floats per branch).
 *
 * @return The starting index in octreeData for this node (branch or leaf).
 */
int buildOctreeRecursive(
    const simd_float3* positions,
    const std::vector<int>& subsetIndices,
    int level,
    int maxLeafSize,
    std::vector<int>& octreeData,
    std::vector<std::vector<int>>& treeLevels,
    int &nodeValues
) {
    // If #particles <= maxLeafSize, create a LEAF node
    if ((int)subsetIndices.size() <= maxLeafSize or level == 5) {
        int nodeStart = (int)octreeData.size();
        // Add this node index to treeLevels
        if (level >= (int)treeLevels.size()) {
            treeLevels.resize(level + 1);
        }
        treeLevels[level].push_back(nodeStart);

        // Write leaf data
        int N = (int)subsetIndices.size();
        octreeData.push_back(N); // number of particles
        octreeData.push_back(nodeValues);
        nodeValues += 1;
        for (int pid : subsetIndices) {
            octreeData.push_back(pid); // each particle index
        }

        return nodeStart;
    }

    // Otherwise, create a BRANCH node
    int nodeStart = (int)octreeData.size();

    // Reserve space for the 10 integers
    for (int i = 0; i < 10; i++) {
        octreeData.push_back(-1);
    }

    // Mark as branch
    octreeData[nodeStart + 0] = 0; // branch marker

    // Store that offset in the second int
    octreeData[nodeStart + 1] = nodeValues;
    nodeValues += 1;

    // Add to levels
    if (level >= (int)treeLevels.size()) {
        treeLevels.resize(level + 1);
    }
    treeLevels[level].push_back(nodeStart);
    
    // Calculate mid point of all particles
    simd_float3 center = { 0, 0, 0 };
    for (int pid : subsetIndices) {
        center += positions[pid];
    }
    center /= subsetIndices.size();

    // Distribute the subset particles into child octants
    std::vector<std::vector<int>> childSubsets(8);
    childSubsets.resize(8);

    for (int pid : subsetIndices) {
        simd_float3 p = positions[pid];
        int octIndex = 0;
        octIndex |= (p.x > center.x) ? 1 : 0;
        octIndex |= (p.y > center.y) ? 2 : 0;
        octIndex |= (p.z > center.z) ? 4 : 0;
        childSubsets[octIndex].push_back(pid);
    }

    // Recursively build each child
    for (int c = 0; c < 8; c++) {
        if (!childSubsets[c].empty()) {
            int childIdx = buildOctreeRecursive(
                positions,
                childSubsets[c],
                level + 1,
                maxLeafSize,
                octreeData,
                treeLevels,
                nodeValues
            );
            // Store child pointer
            octreeData[nodeStart + 2 + c] = childIdx;
        } else {
            // If empty, remains -1
            octreeData[nodeStart + 2 + c] = -1;
        }
    }

    return nodeStart;
}

/**
 * Top-level function to build the octree.
 *
 * @param positions      Pointer to an array of simd_float3 positions.
 * @param positionCount  Number of positions in that array.
 * @param octreeData     (Output) Flat array of integers describing the octree structure.
 * @param treeLevels     (Output) treeLevels[d] = list of node indices at depth d.
 * @param nodeValues     (Output) For each branch node, we store 7 floats:
 *                        [ mass, posX, posY, posZ, gradX, gradY, gradZ ]
 * @param maxLeafSize    Max #particles in a leaf before splitting.
 */
void buildOctree(
    const simd_float3* positions,
    bool* isAlive,
    int positionCount,
    std::vector<int>& octreeData,
    std::vector<std::vector<int>>& treeLevels,
    int& nodeValues,
    int maxLeafSize
) {
    octreeData.clear();
    treeLevels.clear();
    nodeValues = 0;

    if (positionCount <= 0) {
        return;
    }

    // 2. Create a list of all particle indices
    std::vector<int> allIndices;
    for (int i = 0; i < positionCount; i++) {
        if (isAlive[i]) {
            allIndices.push_back(i);
        }
    }

    // 3. Recursively build the octree from the root
    buildOctreeRecursive(
        positions,
        allIndices,
        0,                // level
        maxLeafSize,
        octreeData,
        treeLevels,
        nodeValues
    );
}
