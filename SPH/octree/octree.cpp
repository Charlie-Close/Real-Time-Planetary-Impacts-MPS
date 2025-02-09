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
#include "Parameters.h"

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
    if ((int)subsetIndices.size() <= maxLeafSize or level == MAX_TREE_DEPTH) {
        int nodeStart = (int)octreeData.size();
        // Add this node index to treeLevels
        if (level >= (int)treeLevels.size()) {
            treeLevels.resize(level + 1);
        }
        treeLevels[level].push_back(nodeStart);

        
        // Write leaf data:
        // [number of particles, data pointer (node values), p0, p1, p2, ... pN]
        int N = (int)subsetIndices.size();
        octreeData.push_back(N);
        octreeData.push_back(nodeValues);
        nodeValues += 1;
        for (int pid : subsetIndices) {
            octreeData.push_back(pid);
        }

        return nodeStart;
    }

    // Otherwise, create a BRANCH node
    int nodeStart = (int)octreeData.size();

    // Reserve space for the 10 integers:
    // [0, data pointer (node values), 8xchild pointers]
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

void buildOctree(
    const simd_float3* positions,
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

    std::vector<int> allIndices;
    for (int i = 0; i < positionCount; i++) {
        allIndices.push_back(i);
    }

    buildOctreeRecursive(
        positions,
        allIndices,
        0,
        maxLeafSize,
        octreeData,
        treeLevels,
        nodeValues
    );
}
