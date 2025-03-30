//
//  ANEOS.hpp
//  SPH
//
//  Created by Charlie Close on 27/01/2025.
//

#ifndef ANEOS_hpp
#define ANEOS_hpp

#include <vector>
#include <string>
#include <Metal/Metal.hpp>
#include <simd/simd.h>

struct ANEOSTable
{
    int resolution;

    float* rho;
    float* u;
    float* p;
    float* T;
    
    simd_float4* data; // ( pressure, sound_speed, temperature )
    float* densities;
    
    float minRho;
    float maxRho;
    float minU;
    float maxU;
    float minP;
    float maxP;
    float minT;
    float maxT;
};

ANEOSTable loadANEOSDataFromFile(const std::string &filePath, const int resolution);
ANEOSTable loadHMDataFromFile(const std::string &filePath, const int resolution);
MTL::Texture* createRG32FloatTexture(MTL::Device* device, MTL::CommandQueue* commandQueue, const ANEOSTable& table);

#endif /* ANEOS_hpp */
