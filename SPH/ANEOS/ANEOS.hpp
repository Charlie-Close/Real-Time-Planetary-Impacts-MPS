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
    
    simd_float2* data; // ( pressure, sound_speed )
    
    float minRho;
    float maxRho;
    float minU;
    float maxU;
};

ANEOSTable loadANEOSDataFromFile(const std::string &filePath, const int resolution);
MTL::Texture* createRG32FloatTexture(MTL::Device* device, MTL::CommandQueue* commandQueue, const ANEOSTable& table);

#endif /* ANEOS_hpp */
