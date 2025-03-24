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
#include <thread>

void splitGroupByAxis(const simd_float3* positions, const std::vector<int>& subsetIndices, int axis, std::vector<int>& lower, std::vector<int>& upper) {
    if (subsetIndices.empty()) {
        return;
    }
    
//    simd::float3 min(MAXFLOAT);
//    simd::float3 max(-MAXFLOAT);
//    for (int pid : subsetIndices) {
//        min = simd::min(min, positions[pid]);
//        max = simd::max(max, positions[pid]);
//    }
//    simd::float3 sizes = max - min;
//    axis = 0;
//    if (sizes[1] > sizes[0]) axis = 1;
//    if (sizes[2] > sizes[axis]) axis = 2;
//    float center = 0.5 * (min[axis] + max[axis]);
    
    float min = MAXFLOAT;
    float max = -MAXFLOAT;
    for (int pid : subsetIndices) {
        min = fmin(min, positions[pid][axis]);
        max = fmax(max, positions[pid][axis]);
    }
    float center = 0.5 * (min + max);
    
    long n_lower = 0;
    long n_upper = 0;
    for (int pid : subsetIndices) {
        if (positions[pid][axis] < center) {
            n_lower++;
        } else {
            n_upper++;
        }
    }
    
    lower.reserve(n_lower);
    upper.reserve(n_upper);
    
    for (int pid : subsetIndices) {
        if (positions[pid][axis] < center) {
            lower.push_back(pid);
        } else {
            upper.push_back(pid);
        }
    }
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
    
//    std::vector<int> childSubsets[8];
//    
//    simd::float3 min = simd::float3(MAXFLOAT);
//    simd::float3 max = simd::float3(-MAXFLOAT);
//    for (int pid : subsetIndices) {
//        min = simd::min(min, positions[pid]);
//        max = simd::max(max, positions[pid]);
//    }
//    
//    simd::float3 center = 0.5 * (min + max);
//    
//    long nParticles[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
//    for (int pid : subsetIndices) {
//        uint i = 0;
//        simd::float3 pos = positions[pid];
//        if (pos.x > center.x) i += 1;
//        if (pos.y > center.y) i += 2;
//        if (pos.z > center.z) i += 4;
//        nParticles[i]++;
//    }
//    
//    for (int i = 0; i < 8; i++) {
//        childSubsets[i].reserve(nParticles[i]);
//    }
//    
//    for (int pid : subsetIndices) {
//        uint i = 0;
//        simd::float3 pos = positions[pid];
//        if (pos.x > center.x) i += 1;
//        if (pos.y > center.y) i += 2;
//        if (pos.z > center.z) i += 4;
//        childSubsets[i].push_back(pid);
//    }
    
    std::vector<int> P_1, P_2, P_11, P_12, P_21, P_22, P_111, P_112, P_121, P_211, P_212, P_122, P_221, P_222;
    
    splitGroupByAxis(positions, subsetIndices, 0, P_1, P_2);
    if (subsetIndices.size() < 1e4) {
        splitGroupByAxis(positions, P_1, 1, P_11, P_12);
        splitGroupByAxis(positions, P_2, 1, P_21, P_22);
        splitGroupByAxis(positions, P_11, 2, P_111, P_112);
        splitGroupByAxis(positions, P_12, 2, P_121, P_122);
        splitGroupByAxis(positions, P_21, 2, P_211, P_212);
        splitGroupByAxis(positions, P_22, 2, P_221, P_222);
    } else {
        std::thread t0(splitGroupByAxis, positions, std::ref(P_1), 1, std::ref(P_11), std::ref(P_12));
        std::thread t1(splitGroupByAxis, positions, std::ref(P_2), 1, std::ref(P_21), std::ref(P_22));
        t0.join();
        t1.join();
        std::thread t2(splitGroupByAxis, positions, std::ref(P_11), 2, std::ref(P_111), std::ref(P_112));
        std::thread t3(splitGroupByAxis, positions, std::ref(P_12), 2, std::ref(P_121), std::ref(P_122));
        std::thread t4(splitGroupByAxis, positions, std::ref(P_21), 2, std::ref(P_211), std::ref(P_212));
        std::thread t5(splitGroupByAxis, positions, std::ref(P_22), 2, std::ref(P_221), std::ref(P_222));
        t2.join();
        t3.join();
        t4.join();
        t5.join();
    }
    
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
    int maxLeafSize,
    long estimOctreeDataSize
) {
    octreeData.clear();
    octreeData.reserve(estimOctreeDataSize);
    std::vector<long> treeLevelsSizes;
    for (int i = 0; i < treeLevels.size(); i++) {
        treeLevelsSizes.push_back(treeLevels[i].size());
        treeLevels[i].clear();
    }
    for (int i = 0; i < treeLevels.size(); i++) {
        treeLevels[i].reserve(treeLevelsSizes[i]);
    }
//    treeLevels.clear();
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
    
    for (long i = treeLevels.size() - 1; i >= 0; i--) {
        if (treeLevels[i].size() == 0) {
            treeLevels.erase(treeLevels.begin() + i);
        }
    }
    
    allIndices.clear();
}
