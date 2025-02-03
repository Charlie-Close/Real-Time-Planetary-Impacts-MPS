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

ANEOSTable loadANEOSDataFromFile(const std::string &filePath, const int resolution)
{
    ANEOSTable table;
    table.resolution = resolution;

    std::ifstream infile(filePath);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open ANEOS data file: " + filePath);
    }

    // Skip or parse any header lines if needed
    // e.g. read the version date:
    std::string line;
    for (int i = 0; i < 13; i++) {
        std::getline(infile, line);
    }
    
    int numRho;
    int numT;

    // Next line has num_rho and num_T
    {
        std::getline(infile, line);
        std::istringstream iss(line);
        iss >> numRho >> numT;
    }
    
    float rhoTemp[numRho];
    // Read the rho array
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
    
    float tTemp[numT];
    // Read the T array
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
            throw std::runtime_error("Not enough P data found in file");
        }    }

    infile.close();
    
    table.rho = new float[resolution];
    table.u = new float[resolution];
    table.minRho = 1;
    table.minU = 10000;
    table.maxRho = 29800;
    table.maxU = 8.80519e+10;
    for (int i = 0; i < resolution; i++) {
        table.rho[i] = table.minRho * pow((table.maxRho / table.minRho), (float)i / resolution);
        table.u[i] = table.minU * pow((table.maxU / table.minU), (float)i / resolution);
    }
    
    // RESAMPLE
    table.data = new simd_float2[resolution * resolution];
    int rhoTi = 0;
    for (int rho_i = 0; rho_i < resolution; rho_i++) {
        float rho = table.rho[rho_i];
        while (rho > rhoTemp[rhoTi + 1]) {
            rhoTi++;
        }
        float rhoLFact = (rho - rhoTemp[rhoTi]) / (rhoTemp[rhoTi + 1] - rhoTemp[rhoTi]);
        float rhoHFact = 1 - rhoLFact;
        
        int ulTi = 0;
        int uhTi = 0;
        for (int u_i = 0; u_i < resolution; u_i++) {
            float u = table.u[u_i];
            while (u > uTemp[rhoTi][ulTi + 1]) {
                ulTi++;
            }
            while (u > uTemp[rhoTi + 1][uhTi + 1]) {
                uhTi++;
            }
            
            float uLLFact = (u - uTemp[rhoTi][ulTi]) / (uTemp[rhoTi][ulTi + 1] - uTemp[rhoTi][ulTi]);
            float uLHFact = (u - uTemp[rhoTi + 1][uhTi]) / (uTemp[rhoTi + 1][uhTi + 1] - uTemp[rhoTi + 1][uhTi]);
            float uHLFact = 1 - uLLFact;
            float uHHFact = 1 - uLHFact;
            
            table.data[rho_i * resolution + u_i] = rhoLFact * (dataTemp[rhoTi][ulTi] * uLLFact + dataTemp[rhoTi][ulTi + 1] * uHLFact) +
                                                    rhoHFact * (dataTemp[rhoTi + 1][uhTi] * uLHFact + dataTemp[rhoTi + 1][uhTi + 1] * uHHFact);
        }
    }
    
    return table;
}

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
