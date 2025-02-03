//
//  hdfHandler.hpp
//  SPH
//
//  Created by Charlie Close on 03/11/2024.
//

#ifndef hdfHandler_hpp
#define hdfHandler_hpp

#include <stdio.h>
#include <vector>
#include <simd/simd.h>

struct DataStruct {
    std::vector<simd_float3> positions;
    std::vector<simd_float3> velocities;
    std::vector<float> densities;
    std::vector<float> entropies;
    std::vector<float> internalEnergy;
    std::vector<float> masses;
    std::vector<int> materialIDs;
    std::vector<float> potentials;
    std::vector<float> pressures;
    std::vector<float> smoothingLengths;
};

DataStruct readHDFFile(const std::string& filepath);
void writeHDFFile(const std::string& filepath, const DataStruct& data);

#endif /* hdfHandler_hpp */
