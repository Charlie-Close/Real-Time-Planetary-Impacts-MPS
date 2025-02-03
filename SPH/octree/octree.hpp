//
//  octree.hpp
//  SPH
//
//  Created by Charlie Close on 03/11/2024.
//

#ifndef octree_hpp
#define octree_hpp

#include <vector>
#include <simd/simd.h>

void buildOctree(const simd_float3* positions,
                 bool* isAlive,
                 int positionCount,
                 std::vector<int>& octreeData,
                 std::vector<std::vector<int>>& treeLevels,
                 int& nodeValues,
                 int maxLeafSize = 10);

#endif /* octree_hpp */
