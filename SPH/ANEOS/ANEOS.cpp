//
//  ANEOS.cpp
//  SPH
//
//  Created by Charlie Close on 27/01/2025.
//

#include "ANEOS.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include "Buffers.hpp"
#include "../Parameters.h"

// The data we get in from the ANEOS table is density against temperature to give us internal energy, pressure, sound speed and entropy.
// What we need in our density pass is something which from density and internal energy gives us pressure and sound speed (we don't know
// the temperature). To avoid overly complicated maths in our shaders, we can instead resample our data during initialisation. This function
// reads in data from the file, and produces a table of float2s containing pressure and sound speed, which we can then turn into a texture
// to quickly sample on the GPU.
ANEOSTable loadANEOSDataFromFile(const std::string &filePath, const int resolution)
{
    ANEOSTable table;
    table.resolution = resolution;

    std::ifstream infile(filePath);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open ANEOS data file: " + filePath);
    }

    // Skip header lines
    std::string line;
    for (int i = 0; i < 13; i++) {
        std::getline(infile, line);
    }
    
    // Get the number of densities and temperatures we are going to read.
    int numRho;
    int numT;
    {
        std::getline(infile, line);
        std::istringstream iss(line);
        iss >> numRho >> numT;
    }
    
    // Temporary array for storing the density data (we are going to resample this later)
    float rhoTemp[numRho];
    {
        int count = 0;
        while (count < numRho && std::getline(infile, line)) {
            std::istringstream iss(line);
            float val;
            while (iss >> val && count < numRho) {
                rhoTemp[count++] = val;
                if (count == 0) {
                    table.minRho = val;
                    table.maxRho = val;
                }
                if (val < table.minRho) {
                    table.minRho = val;
                }
                if (val > table.maxRho) {
                    table.maxRho = val;
                }
            }
        }
    }
    
    // Temporary array for storing temperature data (we are going to resample this later)
    float tTemp[numT];
    {
        int count = 0;
        while (count < numT && std::getline(infile, line)) {
            std::istringstream iss(line);
            float val;
            while (iss >> val && count < numT) {
                tTemp[count++] = val;
            }
        }
    }

    // 2D arrays for storing internal energy, and the data we are going to store ( pressure, sound speed )
    float uTemp[numRho][numT];
    simd_float2 dataTemp[numRho][numT];
    {
        int rho_i = 0;
        int t_i = 0;
        while (t_i < numT && std::getline(infile, line)) {
            std::istringstream iss(line);
            float val;
            iss >> val;
            uTemp[rho_i][t_i] = val;
            if (t_i == 0 and rho_i == 0) {
                table.minU = val;
                table.maxU = val;
            }
            if (val < table.minU) {
                table.minU = val;
            }
            if (val > table.maxU) {
                table.maxU = val;
            }
            
            float p;
            float c;
            iss >> p;
            iss >> c;
            dataTemp[rho_i][t_i] = { p * pow(10.f, -18.f), c * pow(10.f, -6.f) };
            
            rho_i++;
            if (rho_i == numRho) {
                rho_i = 0;
                t_i++;
            }
        }
        if (t_i < numT) {
            throw std::runtime_error("Not enough pressure data found in file");
        }
    }
    infile.close();
    
    // Now we can start resampling.
    table.rho = new float[resolution];
    table.u = new float[resolution];
    table.minRho = ANEOS_MIN_RHO * 1e+6;
    table.minU = ANEOS_MIN_U * 1e+12;
    table.maxRho = ANEOS_MAX_RHO * 1e+6;
    table.maxU = ANEOS_MAX_U * 1e+12;
    // We sample exponentially.
    for (int i = 0; i < resolution; i++) {
        table.rho[i] = table.minRho * pow((table.maxRho / table.minRho), (float)i / resolution);
        table.u[i] = table.minU * pow((table.maxU / table.minU), (float)i / resolution);
    }
    
    // We do bilinear interpolation to resample our data
    table.data = new simd_float2[resolution * resolution];
    // Index of which rho we are looking at. Avoids us having to do a binary search each time.
    int rhoTi = 0;
    for (int rho_i = 0; rho_i < resolution; rho_i++) {
        // Find which density indexes we are between
        float rho = table.rho[rho_i];
        while (rho > rhoTemp[rhoTi + 1] and rhoTi < numRho - 2) {
            rhoTi++;
        }
        // Get the interpolation factors (L is lower, H is higher)
        float rhoLFact = (rho - rhoTemp[rhoTi]) / (rhoTemp[rhoTi + 1] - rhoTemp[rhoTi]);
        float rhoHFact = 1 - rhoLFact;
        
        // Now for both the lower and higher indices we need to do the same for internal energy.
        int ulTi = 0;
        int uhTi = 0;
        for (int u_i = 0; u_i < resolution; u_i++) {
            // Find which internal energy indices we are between for both the lower and higher density
            float u = table.u[u_i];
            while (u > uTemp[rhoTi][ulTi + 1] and ulTi < numT - 2) {
                ulTi++;
            }
            while (u > uTemp[rhoTi + 1][uhTi + 1] and ulTi < numT - 2) {
                uhTi++;
            }
            
            // Now we do the bilinear interpolation. Need higher and lower factor for both the higher and lower density.
            float uLLFact = (u - uTemp[rhoTi][ulTi]) / (uTemp[rhoTi][ulTi + 1] - uTemp[rhoTi][ulTi]);
            float uLHFact = (u - uTemp[rhoTi + 1][uhTi]) / (uTemp[rhoTi + 1][uhTi + 1] - uTemp[rhoTi + 1][uhTi]);
            float uHLFact = 1 - uLLFact;
            float uHHFact = 1 - uLHFact;
            
            // Perform the bilinear interpolation and store.
            table.data[rho_i * resolution + u_i] = rhoLFact * (dataTemp[rhoTi][ulTi] * uLLFact + dataTemp[rhoTi][ulTi + 1] * uHLFact) +
                                                    rhoHFact * (dataTemp[rhoTi + 1][uhTi] * uLHFact + dataTemp[rhoTi + 1][uhTi + 1] * uHHFact);
        }
    }
    
    return table;
}

// Turns our ANEOS table into a texture.
MTL::Texture* createRG32FloatTexture(MTL::Device* device, MTL::CommandQueue* commandQueue, const ANEOSTable& table)
{
    using namespace MTL;

    TextureDescriptor* privateDesc = TextureDescriptor::texture2DDescriptor(
        PixelFormatRG32Float,
        table.resolution,
        table.resolution,
        false // mipmapped
    );
    privateDesc->setUsage(TextureUsageShaderRead);
    privateDesc->setStorageMode(StorageModePrivate);

    Texture* privateTexture = device->newTexture(privateDesc);
    
    
    TextureDescriptor* sharedDesc = TextureDescriptor::texture2DDescriptor(
        PixelFormatRG32Float,
        table.resolution,
        table.resolution,
        false // mipmapped
    );
    sharedDesc->setUsage(TextureUsageShaderRead);
    sharedDesc->setStorageMode(StorageModeShared);

    Texture* sharedTexture = device->newTexture(sharedDesc);

    // fill the texture
    Region region = Region::Make2D(0, 0, table.resolution, table.resolution);
    size_t rowBytes = static_cast<size_t>(table.resolution) * sizeof(simd_float2);

    sharedTexture->replaceRegion(region, 0, table.data, rowBytes);
    
    copyDataToTexture(commandQueue, privateTexture, sharedTexture);

    return privateTexture;
}
