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

std::pair<std::vector<int>, std::vector<int>> splitGroupByAxis(const simd_float3* positions, const std::vector<int>& subsetIndices, int axis) {
    std::vector<int> lower, upper;

    if (subsetIndices.empty()) {
        return { lower, upper };
    }
    
    float min = MAXFLOAT;
    float max = - MAXFLOAT;
    for (int pid : subsetIndices) {
        min = fmin(min, positions[pid][axis]);
        max = fmax(max, positions[pid][axis]);
    }
    
    float center = 0.5 * (min + max);
    
    for (int pid : subsetIndices) {
        if (positions[pid][axis] < center) {
            lower.push_back(pid);
        } else {
            upper.push_back(pid);
        }
    }
    
    return { lower, upper };
}

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
    
    // Calculate the size in each dimension:
    simd::float3 min(MAXFLOAT);
    simd::float3 max(-MAXFLOAT);
    for (int pid : subsetIndices) {
        min = simd::min(min, positions[pid]);
        max = simd::max(max, positions[pid]);
    }
    const simd::float3 size = max - min;
    std::array<std::pair<int,float>, 3> axisSize = {{
        {0, size.x},
        {1, size.y},
        {2, size.z}
    }};
    std::sort(axisSize.begin(), axisSize.end(),
              [](auto &a, auto &b){ return a.second > b.second; });
    int axisOrder[3] = {
        axisSize[0].first,
        axisSize[1].first,
        axisSize[2].first
    };
    
    // 4) Split in that order
    auto [P_1, P_2] = splitGroupByAxis(positions, subsetIndices, axisOrder[0]);
    auto [P_11, P_12] = splitGroupByAxis(positions, P_1, axisOrder[1]);
    auto [P_21, P_22] = splitGroupByAxis(positions, P_2, axisOrder[1]);
    auto [P_111, P_112] = splitGroupByAxis(positions, P_11, axisOrder[2]);
    auto [P_121, P_122] = splitGroupByAxis(positions, P_12, axisOrder[2]);
    auto [P_211, P_212] = splitGroupByAxis(positions, P_21, axisOrder[2]);
    auto [P_221, P_222] = splitGroupByAxis(positions, P_22, axisOrder[2]);
    
    std::vector<int> childSubsets[8] = { P_111, P_112, P_121, P_122, P_211, P_212, P_221, P_222 };
    
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
    bool* alive,
    int maxLeafSize
) {
    octreeData.clear();
    for (int i = 0; i < treeLevels.size(); i++) {
        treeLevels[i].clear();
    }
    treeLevels.clear();
    nodeValues = 0;

    if (positionCount <= 0) {
        return;
    }

    std::vector<int> allIndices;
    for (int i = 0; i < positionCount; i++) {
        if (alive[i]) {
            allIndices.push_back(i);
        }
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
    
    allIndices.clear();
}
