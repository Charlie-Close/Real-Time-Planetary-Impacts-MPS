//
//  Compute.hpp
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#ifndef Compute_hpp
#define Compute_hpp

#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <QuartzCore/QuartzCore.hpp>
#include <simd/simd.h>
#include <string>
#include <vector>
#include <AppKit/AppKit.hpp>
#include <MetalKit/MetalKit.hpp>
#include "hdfHandler.hpp"
#include "Parameters.h"

typedef struct {
    simd::float3 pos;
    simd::float3 min;
    simd::float3 max;
    float size;
    float expansion[N_EXPANSION_TERMS];
    float power[P+1];
    float minGrav;
    bool active;
} Multipole;


typedef struct {
    simd::float3 pos;
    float expansion[N_EXPANSION_TERMS];
} Local;


class Compute {
public:
    Compute(MTL::Device* device);
    void updateOctreeBuffer(MTL::Device* device);
    void updateOctreeBufferTest(MTL::Device* device);
    void gravitationalPass(MTL::CommandBuffer* commandBuffer);
    void densityPass(MTL::CommandBuffer* commandBuffer);
    void accelerationPass(MTL::CommandBuffer* commandBuffer);
    void stepPass(MTL::CommandBuffer* commandBuffer);
    
    void sort(MTL::CommandBuffer* commandBuffer);


    // Public buffers (as used by other shaders)
    MTL::Buffer* positionBuffer;
    MTL::Buffer* materialIdBuffer;
    MTL::Buffer* densityBuffer;
    
    int nParticles;

private:
    void buildShaders(MTL::Device* device);
    void buildBuffers(MTL::Device* device, MTL::CommandQueue* commandQueue);
    void loadInitialConditions(MTL::Device* device, MTL::CommandQueue* commandQueue, DataStruct data);
    void updateOctreeData(MTL::Device* device);
    void encodeCommand(MTL::ComputeCommandEncoder* computeEncoder, MTL::ComputePipelineState* command, long size);
    
    std::vector<std::vector<int>> treeLevels;
    std::vector<std::vector<int>> treeLevelsTemp;
    std::vector<int> octreeData;
    int nodeValues;
    bool updating = false;
    int nBlocks;
    long prevGravDataSize = 0;
    long prevTreeValuesSize = 0;
    
    int cellsPerDim;
    float cellSize;

    MTL::ComputePipelineState* _densityPSO;
    MTL::ComputePipelineState* _accelerationPSO;
    MTL::ComputePipelineState* _gravityPSO;
    MTL::ComputePipelineState* _upTreePSO;
    MTL::ComputePipelineState* _downTreePSO;
    MTL::ComputePipelineState* _mStepPSO;
    
    MTL::ComputePipelineState* _hashPSO;
    MTL::ComputePipelineState* _histPSO;
    MTL::ComputePipelineState* _scanPSO;
    MTL::ComputePipelineState* _sumPSO;
    MTL::ComputePipelineState* _sortPSO;
    MTL::ComputePipelineState* _initialisePSO;
    MTL::ComputePipelineState* _findPosPSO;


    
    MTL::Buffer* _velocityBuffer;
    MTL::Buffer* _accelerationBuffer;
    MTL::Buffer* _internalEnergyBuffer;
    MTL::Buffer* _massBuffer;
    MTL::Buffer* _pressureBuffer;
    MTL::Buffer* _smoothingLengthBuffer;
    MTL::Buffer* _gradientTermsBuffer;
    MTL::Buffer* _speedOfSoundBuffer;
    MTL::Buffer* _dInternalEnergyBuffer;
    MTL::Buffer* _potentialGradientBuffer;
    MTL::Buffer* _cellArrayi;
    MTL::Buffer* _cellArrayj;
    MTL::Buffer* _bucketHist;
    MTL::Buffer* _bucketOffset;
    MTL::Buffer* _particleOffset;
    MTL::Buffer* _ittr[SORTING_ITTERATIONS];
    MTL::Buffer* _nParticles;
    MTL::Buffer* _nBlocks;
    MTL::Buffer* _cellStart;
    MTL::Buffer* _cellEnd;
    MTL::Buffer* _cellSize;
    MTL::Buffer* _cellsPerDim;
    MTL::Buffer* _tree;
    MTL::Buffer* _treeTmp;
    MTL::Buffer* _multipoleExpansions;
    MTL::Buffer* _localExpansion;
    MTL::Buffer* _parentIndexes;
    MTL::Buffer* _leafPointers;
    MTL::Buffer* _localGravi;
    MTL::Buffer* _localGravj;
    MTL::Buffer* _active;
    MTL::Buffer* _gravAbs;
    std::vector<MTL::Buffer*> _treeLevelBuffers;
    MTL::Buffer* _dt;
    MTL::Buffer* _balsara;
    MTL::Buffer* _nextActiveTime;
    MTL::Buffer* _globalTime;
    MTL::Buffer* _dhdt;

    
    MTL::Texture* _forsterite;
    MTL::Texture* _Fe85Si15;
};

#endif /* Compute_hpp */
