//
//  Compute.cpp
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include "Compute.hpp"
#include "Buffers.hpp"
#include <iostream>
#include "octree.hpp"
#include "ANEOS.hpp"
#include <thread>
#include "lodepng.h"

// -------------------------------- //
//                                  //
//          Constructor             //
//                                  //
// -------------------------------- //

Compute::Compute(MTL::Device* device) {
    DataStruct data = readHDFFile(FILEPATH);
    
    MTL::CommandQueue* commandQueue = device->newCommandQueue();
        
    nParticles = (int)data.positions.size();
    cellSize = CELL_WIDTH;
    cellsPerDim = pow(2, CELL_POWER);
    
    std::cout << nParticles << std::endl;

    buildShaders(device);
    buildBuffers(device, commandQueue);
    loadInitialConditions(device, commandQueue, data);
    updateOctreeData(device);
//    updateOctreeBuffer(device);
    
    // Sort data
    MTL::CommandBuffer* cmdBuffer = commandQueue->commandBuffer();
    sort(cmdBuffer);
    cmdBuffer->commit();
    cmdBuffer->waitUntilCompleted();
    shuffleData(device, commandQueue);
}

Compute::~Compute() {
    positionBuffer->release();
    materialIdBuffer->release();
    densityBuffer->release();
    _velocityBuffer->release();
    _accelerationBuffer->release();
    _internalEnergyBuffer->release();
    _massBuffer->release();
    _pressureBuffer->release();
    _smoothingLengthBuffer->release();
    _gradientTermsBuffer->release();
    _speedOfSoundBuffer->release();
    _dInternalEnergyBuffer->release();
    _cellArrayi->release();
    _cellArrayj->release();
    _bucketHist->release();
    _bucketOffset->release();
    _particleOffset->release();
    for (int i = 0; i < SORTING_ITTERATIONS; i++) {
        _ittr[i]->release();
    }
    _nParticles->release();
    _nBlocks->release();
    _cellSize->release();
    _cellEnd->release();
    _cellSize->release();
    _cellsPerDim->release();
    _tree->release();
    _multipoleExpansions->release();
    _localExpansion->release();
    _parentIndexes->release();
    _leafPointers->release();
    _localGravi->release();
    _localGravj->release();
    _active->release();
    _gravAbs->release();
    for (int i = 0; i < _treeLevelBuffers.size(); i++) {
        _treeLevelBuffers[i]->release();
    }
    _dt->release();
    _balsara->release();
    _nextActiveTime->release();
    _globalTime->release();
    _dhdt->release();
    _forsterite->release();
    _Fe85Si15->release();
    
    _hashPSO->release();
    _histPSO->release();
    _scanPSO->release();
    _sumPSO->release();
    _sortPSO->release();
    _initialisePSO->release();
    _findPosPSO->release();
    _upTreePSO->release();
    _downTreePSO->release();
    _densityPSO->release();
    _accelerationPSO->release();
    _mStepPSO->release();
    _shuffleFloat->release();
    _shuffleFloat3->release();
    _shuffleInt->release();
    _shuffleTree->release();
}

// -------------------------------- //
//                                  //
//   Initial Condition Generation   //
//                                  //
// -------------------------------- //

void Compute::loadInitialConditions(MTL::Device* device, MTL::CommandQueue* commandQueue, DataStruct data) {
    ANEOSTable forTable = loadANEOSDataFromFile("ANEOS_forsterite_S19.txt", ANEOS_TEXTURE_RESOLUTION);
    ANEOSTable feTable = loadANEOSDataFromFile("ANEOS_Fe85Si15_S20.txt", ANEOS_TEXTURE_RESOLUTION);
    
    _forsterite = createRG32FloatTexture(device, commandQueue, forTable);
    _Fe85Si15 = createRG32FloatTexture(device, commandQueue, feTable);
    
    // Write all the data to the buffers
    writeDataToBuffer(positionBuffer, data.positions);
    writeDataToPrivateBuffer(device, commandQueue, _velocityBuffer, data.velocities);
    writeDataToPrivateBuffer(device, commandQueue, densityBuffer, data.densities);
    writeDataToPrivateBuffer(device, commandQueue, _internalEnergyBuffer, data.internalEnergy);
    writeDataToPrivateBuffer(device, commandQueue, _massBuffer, data.masses);
    writeDataToPrivateBuffer(device, commandQueue, _pressureBuffer, data.pressures);
    writeDataToPrivateBuffer(device, commandQueue, _smoothingLengthBuffer, data.smoothingLengths);
    writeDataToPrivateBuffer(device, commandQueue, materialIdBuffer, data.materialIDs);
    writeDataToPrivateBuffer(device, commandQueue, _nParticles, nParticles);
    for (uint i = 0; i < SORTING_ITTERATIONS; i++) {
        writeDataToPrivateBuffer(device, commandQueue, _ittr[i], i);
    }
    writeDataToPrivateBuffer(device, commandQueue, _nBlocks, nBlocks);
    
    bool* active = new bool[nParticles];
    bool* alive = new bool[nParticles];
    for (int i = 0; i < nParticles; i++) {
        active[i] = true;
        alive[i] = true;
    }
    writeDataToPrivateBuffer(device, commandQueue, _active, active, nParticles);
    writeDataToPrivateBuffer(device, commandQueue, _alive, alive, nParticles);
    writeDataToPrivateBuffer(device, commandQueue, _globalTime, (int)0);
    
    int* nextActiveTime = new int[nParticles];
    for (int i = 0; i < nParticles; i++) {
        nextActiveTime[i] = 0;
    }
    writeDataToPrivateBuffer(device, commandQueue, _nextActiveTime, nextActiveTime, nParticles);
    
    writeDataToPrivateBuffer(device, commandQueue, _cellSize, cellSize);
    writeDataToPrivateBuffer(device, commandQueue, _cellsPerDim, cellsPerDim);
}

// -------------------------------- //
//                                  //
//      Metal Object Builders       //
//                                  //
// -------------------------------- //

void Compute::buildBuffers(MTL::Device* device, MTL::CommandQueue *commandQueue) {
    nBlocks = (nParticles / SORTING_BLOCK_SIZE) + 1;
    std::cout << "N BLOCK " << nBlocks << std::endl;
    positionBuffer = device->newBuffer(nParticles * sizeof(simd_float3), MTL::ResourceStorageModeShared);
    _velocityBuffer = device->newBuffer(nParticles * sizeof(simd_float3), MTL::ResourceStorageModePrivate);
    _accelerationBuffer = device->newBuffer(nParticles * sizeof(simd_float3), MTL::ResourceStorageModePrivate);
    densityBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _internalEnergyBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _massBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _pressureBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _smoothingLengthBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    materialIdBuffer = device->newBuffer(nParticles * sizeof(int), MTL::ResourceStorageModePrivate);
    _gradientTermsBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _speedOfSoundBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _dInternalEnergyBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _balsara = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    rhoGrads = device->newBuffer(nParticles * sizeof(simd_float3), MTL::ResourceStorageModePrivate);
    _cellArrayi = device->newBuffer(nParticles * sizeof(simd_uint2), MTL::ResourceStorageModePrivate);
    _cellArrayj = device->newBuffer(nParticles * sizeof(simd_uint2), MTL::ResourceStorageModePrivate);
    _bucketHist = device->newBuffer(nBlocks * SORTING_BUCKET_NUMBER * sizeof(uint), MTL::ResourceStorageModePrivate);
    _bucketOffset = device->newBuffer(SORTING_BUCKET_NUMBER * sizeof(uint), MTL::ResourceStorageModePrivate);
    _particleOffset = device->newBuffer(nParticles * sizeof(uint), MTL::ResourceStorageModePrivate);
    for (int i = 0; i < SORTING_ITTERATIONS; i++) {
        _ittr[i] = device->newBuffer(sizeof(uint), MTL::ResourceStorageModePrivate);
    }
    _nParticles = device->newBuffer(sizeof(int), MTL::ResourceStorageModePrivate);
    _nBlocks = device->newBuffer(sizeof(int), MTL::ResourceStorageModePrivate);
    _cellStart = device->newBuffer(cellsPerDim * cellsPerDim * cellsPerDim * sizeof(int), MTL::ResourceStorageModePrivate);
    _cellEnd = device->newBuffer(cellsPerDim * cellsPerDim * cellsPerDim * sizeof(int), MTL::ResourceStorageModePrivate);
    _cellSize = device->newBuffer(sizeof(float), MTL::ResourceStorageModePrivate);
    _cellsPerDim = device->newBuffer(sizeof(int), MTL::ResourceStorageModePrivate);
    _leafPointers = device->newBuffer(nParticles * sizeof(int), MTL::ResourceStorageModePrivate);
    _active = device->newBuffer(sizeof(bool) * nParticles, MTL::ResourceStorageModePrivate);
    _alive = device->newBuffer(sizeof(bool) * nParticles, MTL::ResourceStorageModeShared);
    _dt = device->newBuffer(sizeof(float), MTL::ResourceStorageModePrivate);
    _gravAbs = device->newBuffer(sizeof(float) * nParticles, MTL::ResourceStorageModePrivate);
    _tree = nullptr;
    _nextActiveTime = device->newBuffer(sizeof(uint) * nParticles, MTL::ResourceStorageModePrivate);
    _globalTime = device->newBuffer(sizeof(int), MTL::ResourceStorageModeShared);
    _dhdt = device->newBuffer(sizeof(float) * nParticles, MTL::ResourceStorageModePrivate);
}

void Compute::buildShaders(MTL::Device* device) {
    NS::Error** error = nil;
    
    // Load the shader files with a .metal file extension in the project
    MTL::Library* defaultLibrary = device->newDefaultLibrary();
    
    MTL::Function* densityFunction = defaultLibrary->newFunction(NS::String::string("density", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* accelerationFunction = defaultLibrary->newFunction(NS::String::string("acceleration", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* upPassFunction = defaultLibrary->newFunction(NS::String::string("upPass", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* downPassFunction = defaultLibrary->newFunction(NS::String::string("downPass", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* stepFunction = defaultLibrary->newFunction(NS::String::string("step", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* hashFunction = defaultLibrary->newFunction(NS::String::string("hash", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* histFunction = defaultLibrary->newFunction(NS::String::string("hist", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* scanFunction = defaultLibrary->newFunction(NS::String::string("scan", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* sumFunction = defaultLibrary->newFunction(NS::String::string("sum", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* sortFunction = defaultLibrary->newFunction(NS::String::string("sort", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* initFunction = defaultLibrary->newFunction(NS::String::string("initialise", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* fcpFunction = defaultLibrary->newFunction(NS::String::string("findCellPositions", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* shuffleFloat = defaultLibrary->newFunction(NS::String::string("shuffleFloat", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* shuffleFloat3 = defaultLibrary->newFunction(NS::String::string("shuffleFloat3", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* shuffleInt = defaultLibrary->newFunction(NS::String::string("shuffleInt", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* shuffleBool = defaultLibrary->newFunction(NS::String::string("shuffleBool", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* shuffleTree = defaultLibrary->newFunction(NS::String::string("shuffleTree", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* inverseArray = defaultLibrary->newFunction(NS::String::string("inverseArray", NS::StringEncoding::UTF8StringEncoding));
    
    
    
    // Create a compute pipeline state object.
    _densityPSO = device->newComputePipelineState(densityFunction, error);
    _accelerationPSO = device->newComputePipelineState(accelerationFunction, error);
    _upTreePSO = device->newComputePipelineState(upPassFunction, error);
    _downTreePSO = device->newComputePipelineState(downPassFunction, error);
    _mStepPSO = device->newComputePipelineState(stepFunction, error);
    _hashPSO = device->newComputePipelineState(hashFunction, error);
    _histPSO = device->newComputePipelineState(histFunction, error);
    _scanPSO = device->newComputePipelineState(scanFunction, error);
    _sumPSO = device->newComputePipelineState(sumFunction, error);
    _sortPSO = device->newComputePipelineState(sortFunction, error);
    _initialisePSO = device->newComputePipelineState(initFunction, error);
    _findPosPSO = device->newComputePipelineState(fcpFunction, error);
    _shuffleFloat = device->newComputePipelineState(shuffleFloat, error);
    _shuffleFloat3 = device->newComputePipelineState(shuffleFloat3, error);
    _shuffleInt = device->newComputePipelineState(shuffleInt, error);
    _shuffleBool = device->newComputePipelineState(shuffleBool, error);
    _shuffleTree = device->newComputePipelineState(shuffleTree, error);
    _inverseArray = device->newComputePipelineState(inverseArray, error);
}

// -------------------------------- //
//                                  //
//          Compute Command         //
//                                  //
// -------------------------------- //

void Compute::densityPass(MTL::CommandBuffer* commandBuffer) {
    MTL::ComputeCommandEncoder* computeEncoder = commandBuffer->computeCommandEncoder();
    assert(computeEncoder != nil);
    
    computeEncoder->setBuffer(positionBuffer, 0, 0);
    computeEncoder->setBuffer(_velocityBuffer, 0, 1);
    computeEncoder->setBuffer(densityBuffer, 0, 2);
    computeEncoder->setBuffer(_internalEnergyBuffer, 0, 3);
    computeEncoder->setBuffer(_massBuffer, 0, 4);
    computeEncoder->setBuffer(_pressureBuffer, 0, 5);
    computeEncoder->setBuffer(_smoothingLengthBuffer, 0, 6);
    computeEncoder->setBuffer(materialIdBuffer, 0, 7);
    computeEncoder->setBuffer(_gradientTermsBuffer, 0, 8);
    computeEncoder->setBuffer(_speedOfSoundBuffer, 0, 9);
    computeEncoder->setBuffer(_balsara, 0, 10);
    computeEncoder->setBuffer(rhoGrads, 0, 11);
    computeEncoder->setBuffer(_cellArrayi, 0, 12);
    computeEncoder->setBuffer(_cellStart, 0, 13);
    computeEncoder->setBuffer(_cellEnd, 0, 14);
    computeEncoder->setBuffer(_cellSize, 0, 15);
    computeEncoder->setBuffer(_cellsPerDim, 0, 16);
    computeEncoder->setBuffer(_active, 0, 17);
    computeEncoder->setBuffer(_alive, 0, 18);
    computeEncoder->setBuffer(_dt, 0, 19);
    computeEncoder->setTexture(_forsterite, 0);
    computeEncoder->setTexture(_Fe85Si15, 1);
    
    encodeCommand(computeEncoder, _densityPSO, nParticles, 256);

    //  End the compute pass.
    computeEncoder->endEncoding();
}

void Compute::accelerationPass(MTL::CommandBuffer* commandBuffer) {
    MTL::ComputeCommandEncoder* computeEncoder = commandBuffer->computeCommandEncoder();
    assert(computeEncoder != nil);
    
    computeEncoder->setBuffer(positionBuffer, 0, 0);
    computeEncoder->setBuffer(_velocityBuffer, 0, 1);
    computeEncoder->setBuffer(_accelerationBuffer, 0, 2);
    computeEncoder->setBuffer(densityBuffer, 0, 3);
    computeEncoder->setBuffer(_internalEnergyBuffer, 0, 4);
    computeEncoder->setBuffer(_massBuffer, 0, 5);
    computeEncoder->setBuffer(_pressureBuffer, 0, 6);
    computeEncoder->setBuffer(_smoothingLengthBuffer, 0, 7);
    computeEncoder->setBuffer(_gradientTermsBuffer, 0, 8);
    computeEncoder->setBuffer(_speedOfSoundBuffer, 0, 9);
    computeEncoder->setBuffer(_dInternalEnergyBuffer, 0, 10);
    computeEncoder->setBuffer(_balsara, 0, 11);
    computeEncoder->setBuffer(_dhdt, 0, 12);
    computeEncoder->setBuffer(_cellArrayi, 0, 13);
    computeEncoder->setBuffer(_cellStart, 0, 14);
    computeEncoder->setBuffer(_cellEnd, 0, 15);
    computeEncoder->setBuffer(_cellSize, 0, 16);
    computeEncoder->setBuffer(_cellsPerDim, 0, 17);
    computeEncoder->setBuffer(_active, 0, 18);
    computeEncoder->setBuffer(_alive, 0, 19);
    computeEncoder->setBuffer(_nextActiveTime, 0, 20);
    computeEncoder->setBuffer(_globalTime, 0, 21);
    computeEncoder->setBuffer(_dt, 0, 22);
    encodeCommand(computeEncoder, _accelerationPSO, nParticles, 256);

    //  End the compute pass.
    computeEncoder->endEncoding();
}

void Compute::stepPass(MTL::CommandBuffer* commandBuffer) {
    // Start a compute pass.
    MTL::ComputeCommandEncoder* computeEncoder = commandBuffer->computeCommandEncoder();
    assert(computeEncoder != nil);
    
    computeEncoder->setBuffer(positionBuffer, 0, 0);
    computeEncoder->setBuffer(_velocityBuffer, 0, 1);
    computeEncoder->setBuffer(_accelerationBuffer, 0, 2);
    computeEncoder->setBuffer(densityBuffer, 0, 3);
    computeEncoder->setBuffer(_internalEnergyBuffer, 0, 4);
    computeEncoder->setBuffer(_dInternalEnergyBuffer, 0, 5);
    computeEncoder->setBuffer(_dhdt, 0, 6);
    computeEncoder->setBuffer(_smoothingLengthBuffer, 0, 7);
    computeEncoder->setBuffer(_active, 0, 8);
    computeEncoder->setBuffer(_nextActiveTime, 0, 9);
    computeEncoder->setBuffer(_globalTime, 0, 10);
    computeEncoder->setBuffer(_dt, 0, 11);
    computeEncoder->setBuffer(_alive, 0, 12);
    encodeCommand(computeEncoder, _mStepPSO, nParticles);

    //  End the compute pass.
    computeEncoder->endEncoding();
}

void Compute::sort(MTL::CommandBuffer *commandBuffer) {
    framesSinceLastShuffled++;
    // Calcualte particle indexes and write to cell data
    MTL::ComputeCommandEncoder* computeEncoder = commandBuffer->computeCommandEncoder();
    computeEncoder->setBuffer(positionBuffer, 0, 0);
    computeEncoder->setBuffer(_cellSize, 0, 1);
    computeEncoder->setBuffer(_cellsPerDim, 0, 2);
    computeEncoder->setBuffer(_cellArrayi, 0, 3);
    encodeCommand(computeEncoder, _hashPSO, nParticles);
    computeEncoder->endEncoding();
    
    for (int i = 0; i < SORTING_ITTERATIONS; i++) {
        
        MTL::ComputeCommandEncoder* cEncHist = commandBuffer->computeCommandEncoder();
        cEncHist->setBuffer((i % 2 == 0 ? _cellArrayi : _cellArrayj), 0, 0);
        cEncHist->setBuffer(_bucketHist, 0, 1);
        cEncHist->setBuffer(_particleOffset, 0, 2);
        cEncHist->setBuffer(_ittr[i], 0, 3);
        cEncHist->setBuffer(_nParticles, 0, 4);
        encodeCommand(cEncHist, _histPSO, nBlocks);
        cEncHist->endEncoding();
        
        MTL::ComputeCommandEncoder* cEncScan = commandBuffer->computeCommandEncoder();
        cEncScan->setBuffer((i % 2 == 0 ? _cellArrayi : _cellArrayj), 0, 0);
        cEncScan->setBuffer(_bucketHist, 0, 1);
        cEncScan->setBuffer(_bucketOffset, 0, 2);
        cEncScan->setBuffer(_particleOffset, 0, 3);
        cEncScan->setBuffer(_ittr[i], 0, 4);
        cEncScan->setBuffer(_nBlocks, 0, 5);
        encodeCommand(cEncScan, _scanPSO, SORTING_BUCKET_NUMBER);
        cEncScan->endEncoding();
        
        MTL::ComputeCommandEncoder* cEncSum = commandBuffer->computeCommandEncoder();
        cEncSum->setBuffer(_bucketHist, 0, 0);
        cEncSum->setBuffer(_bucketOffset, 0, 1);
        encodeCommand(cEncSum, _sumPSO, 1);
        cEncSum->endEncoding();
        
        MTL::ComputeCommandEncoder* cEncSort = commandBuffer->computeCommandEncoder();
        cEncSort->setBuffer((i % 2 == 0 ? _cellArrayi : _cellArrayj), 0, 0);
        cEncSort->setBuffer((i % 2 == 0 ? _cellArrayj : _cellArrayi), 0, 1);
        cEncSort->setBuffer(_bucketHist, 0, 2);
        cEncSort->setBuffer(_bucketOffset, 0, 3);
        cEncSort->setBuffer(_particleOffset, 0, 4);
        cEncSort->setBuffer(_ittr[i], 0, 5);
        encodeCommand(cEncSort, _sortPSO, nParticles);
        cEncSort->endEncoding();
    }
    
    MTL::ComputeCommandEncoder* cEncInit = commandBuffer->computeCommandEncoder();
    cEncInit->setBuffer(_cellStart, 0, 0);
    cEncInit->setBuffer(_cellEnd, 0, 1);
    encodeCommand(cEncInit, _initialisePSO, cellsPerDim * cellsPerDim * cellsPerDim);
    cEncInit->endEncoding();

    MTL::ComputeCommandEncoder* cEncPos = commandBuffer->computeCommandEncoder();
    cEncPos->setBuffer(_cellArrayi, 0, 0);
    cEncPos->setBuffer(_cellStart, 0, 1);
    cEncPos->setBuffer(_cellEnd, 0, 2);
    encodeCommand(cEncPos, _findPosPSO, nParticles);
    cEncPos->endEncoding();
}

void Compute::gravitationalPass(MTL::CommandBuffer* commandBuffer) {
    // Go up the tree first:
    for (long i = _treeLevelBuffers.size() - 1; i >= 0; i--) {
        MTL::ComputeCommandEncoder* computeEncoder = commandBuffer->computeCommandEncoder();
        assert(computeEncoder != nil);
        
        computeEncoder->setBuffer(positionBuffer, 0, 0);
        computeEncoder->setBuffer(_massBuffer, 0, 1);
        computeEncoder->setBuffer(_smoothingLengthBuffer, 0, 2);
        computeEncoder->setBuffer(_tree, 0, 3);
        computeEncoder->setBuffer(_multipoleExpansions, 0, 4);
        computeEncoder->setBuffer(_localExpansion, 0, 5);
        computeEncoder->setBuffer(_treeLevelBuffers[i], 0, 6);
        computeEncoder->setBuffer(_parentIndexes, 0, 7);
        computeEncoder->setBuffer(_gravAbs, 0, 8);
        computeEncoder->setBuffer(_active, 0, 9);
        computeEncoder->setBuffer(_nextActiveTime, 0, 10);
        computeEncoder->setBuffer(_globalTime, 0, 11);


        encodeCommand(computeEncoder, _upTreePSO, treeLevels[i].size());

        //  End the compute pass.
        computeEncoder->endEncoding();
    }
    
    for (int i = 0; i < _treeLevelBuffers.size(); i++) {
        MTL::ComputeCommandEncoder* computeEncoder = commandBuffer->computeCommandEncoder();
        assert(computeEncoder != nil);
        
    
        computeEncoder->setBuffer(positionBuffer, 0, 0);
        computeEncoder->setBuffer(_accelerationBuffer, 0, 1);
        computeEncoder->setBuffer(_massBuffer, 0, 2);
        computeEncoder->setBuffer(_smoothingLengthBuffer, 0, 3);
        computeEncoder->setBuffer(_tree, 0, 4);
        computeEncoder->setBuffer(_multipoleExpansions, 0, 5);
        computeEncoder->setBuffer(_localExpansion, 0, 6);
        computeEncoder->setBuffer(_treeLevelBuffers[i], 0, 7);
        computeEncoder->setBuffer(_parentIndexes, 0, 8);
        computeEncoder->setBuffer(i % 2 == 0 ? _localGravi : _localGravj, 0, 9);
        computeEncoder->setBuffer(i % 2 == 0 ? _localGravj : _localGravi, 0, 10);
        computeEncoder->setBuffer(_gravAbs, 0, 11);
        computeEncoder->setBuffer(_active, 0, 12);


        encodeCommand(computeEncoder, _downTreePSO, treeLevels[i].size());

        //  End the compute pass.
        computeEncoder->endEncoding();
    }
}

void Compute::encodeCommand(MTL::ComputeCommandEncoder* computeEncoder, MTL::ComputePipelineState* command, long size, int cap) {
    MTL::Size gridSize = MTL::Size(size, 1, 1);

    computeEncoder->setComputePipelineState(command);
    // Calculate a threadgroup size.
    NS::UInteger threadGroupSize = cap == 0 ? command->maxTotalThreadsPerThreadgroup() : cap;
    if (threadGroupSize > size)
    {
        threadGroupSize = size;
    }
    MTL::Size threadgroupSize = MTL::Size(threadGroupSize, 1, 1);

    // Encode the compute command.
    computeEncoder->dispatchThreads(gridSize, threadgroupSize);
}

float Compute::getTime() {
    int* globalTime = static_cast<int*>(_globalTime->contents());
    return (*globalTime) * DT_MIN1;
}

// -------------------------------- //
//                                  //
//        Data Organisation         //
//                                  //
// -------------------------------- //

void Compute::updateOctreeData(MTL::Device* device) {
    simd_float3* positions = static_cast<simd_float3*>(positionBuffer->contents());
    bool* alive = static_cast<bool*>(_alive->contents());
    buildOctree(positions, nParticles, octreeData, treeLevelsTemp, nodeValues, alive, MAX_CHILDREN_IN_LEAF);
    _treeTmp = device->newBuffer(octreeData.size() * sizeof(int), MTL::ResourceStorageModeShared);
    writeDataToBuffer(_treeTmp, octreeData);
    updating = false;
}

void Compute::organisation(MTL::Device* device, MTL::CommandQueue* commandQueue) {
    if (framesSinceLastShuffled > SHUFFLE_FRAMES and !updating) {
        shuffleData(device, commandQueue);
        framesSinceLastShuffled = 0;
        return;
    }
    if (!updating) {
        updateOctreeBuffer(device);
    }
    // Otherwise nothing to do.
}


void Compute::shuffleData(MTL::Device* device, MTL::CommandQueue* commandQueue) {
    positionBuffer = shuffleFloat3(device, commandQueue, positionBuffer);
    _velocityBuffer = shuffleFloat3(device, commandQueue, _velocityBuffer);
    _accelerationBuffer = shuffleFloat3(device, commandQueue, _accelerationBuffer);
    densityBuffer = shuffleFloat(device, commandQueue, densityBuffer);
    _internalEnergyBuffer = shuffleFloat(device, commandQueue, _internalEnergyBuffer);
    _massBuffer = shuffleFloat(device, commandQueue, _massBuffer);
    _pressureBuffer = shuffleFloat(device, commandQueue, _pressureBuffer);
    _smoothingLengthBuffer = shuffleFloat(device, commandQueue, _smoothingLengthBuffer);
    materialIdBuffer = shuffleInt(device, commandQueue, materialIdBuffer);
    _gradientTermsBuffer = shuffleFloat(device, commandQueue, _gradientTermsBuffer);
    _speedOfSoundBuffer = shuffleFloat(device, commandQueue, _speedOfSoundBuffer);
    _dInternalEnergyBuffer = shuffleFloat(device, commandQueue, _dInternalEnergyBuffer);
    _balsara = shuffleFloat(device, commandQueue, _balsara);
    rhoGrads = shuffleFloat3(device, commandQueue, rhoGrads);
    _gravAbs = shuffleFloat(device, commandQueue, _gravAbs);
    _nextActiveTime = shuffleInt(device, commandQueue, _nextActiveTime);
    _dhdt = shuffleFloat(device, commandQueue, _dhdt);
    _alive = shuffleBool(device, commandQueue, _alive);
    updateOctreeBuffer(device);
    inverseCellParticles(device, commandQueue);
    for (int i = 0; i < _treeLevelBuffers.size(); i++) {
        shuffleTreeLayer(device, commandQueue, i);
    }
}

MTL::Buffer* Compute::shuffleFloat(MTL::Device* device, MTL::CommandQueue* commandQueue, MTL::Buffer* buffer) {
    if (buffer->length() != nParticles * sizeof(float)) {
        std::cout << "uh oh" << std::endl;
    }
    MTL::Buffer* sortedBuffer = device->newBuffer(nParticles * sizeof(float), buffer->resourceOptions());
    bool success = false;
    MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
    cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);
    while (not success) {
        MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer(cmdDesc);
        MTL::ComputeCommandEncoder* computeCommandEncoder = commandBuffer->computeCommandEncoder();
        computeCommandEncoder->setBuffer(buffer, 0, 0);
        computeCommandEncoder->setBuffer(sortedBuffer, 0, 1);
        computeCommandEncoder->setBuffer(_cellArrayi, 0, 2);
        encodeCommand(computeCommandEncoder, _shuffleFloat, nParticles);
        computeCommandEncoder->endEncoding();
        commandBuffer->commit();
        commandBuffer->waitUntilCompleted();
        success = commandBuffer->status() == MTL::CommandBufferStatusCompleted;
        if (!success) {
            std::cout << "Error, retrying" << std::endl;
        }
    }
    cmdDesc->release();
    buffer->release();
    return sortedBuffer;
}

MTL::Buffer* Compute::shuffleFloat3(MTL::Device* device, MTL::CommandQueue* commandQueue, MTL::Buffer* buffer) {
    if (buffer->length() != nParticles * sizeof(simd_float3)) {
        std::cout << "uh oh" << std::endl;
    }
    MTL::Buffer* sortedBuffer = device->newBuffer(nParticles * sizeof(simd_float3), buffer->resourceOptions());
    bool success = false;
    MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
    cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);

    while (not success) {
        MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer(cmdDesc);
        MTL::ComputeCommandEncoder* computeCommandEncoder = commandBuffer->computeCommandEncoder();
        computeCommandEncoder->setBuffer(buffer, 0, 0);
        computeCommandEncoder->setBuffer(sortedBuffer, 0, 1);
        computeCommandEncoder->setBuffer(_cellArrayi, 0, 2);
        encodeCommand(computeCommandEncoder, _shuffleFloat3, nParticles);
        computeCommandEncoder->endEncoding();
        commandBuffer->commit();
        commandBuffer->waitUntilCompleted();
        success = commandBuffer->status() == MTL::CommandBufferStatusCompleted;
        if (!success) {
            std::cout << "Error, retrying" << std::endl;
        }
    }
    cmdDesc->release();
    buffer->release();
    return sortedBuffer;
}

MTL::Buffer* Compute::shuffleInt(MTL::Device* device, MTL::CommandQueue* commandQueue, MTL::Buffer* buffer) {
    if (buffer->length() != nParticles * sizeof(int)) {
        std::cout << "uh oh" << std::endl;
    }
    MTL::Buffer* sortedBuffer = device->newBuffer(nParticles * sizeof(int), buffer->resourceOptions());
    bool success = false;
    MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
    cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);

    while (not success) {
        MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer(cmdDesc);
        MTL::ComputeCommandEncoder* computeCommandEncoder = commandBuffer->computeCommandEncoder();
        computeCommandEncoder->setBuffer(buffer, 0, 0);
        computeCommandEncoder->setBuffer(sortedBuffer, 0, 1);
        computeCommandEncoder->setBuffer(_cellArrayi, 0, 2);
        encodeCommand(computeCommandEncoder, _shuffleInt, nParticles);
        computeCommandEncoder->endEncoding();
        commandBuffer->commit();
        commandBuffer->waitUntilCompleted();
        success = commandBuffer->status() == MTL::CommandBufferStatusCompleted;
        if (!success) {
            std::cout << "Error, retrying" << std::endl;
        }
    }
    cmdDesc->release();
    buffer->release();
    return sortedBuffer;
}

MTL::Buffer* Compute::shuffleBool(MTL::Device* device, MTL::CommandQueue* commandQueue, MTL::Buffer* buffer) {
    if (buffer->length() != nParticles * sizeof(bool)) {
        std::cout << "uh oh" << std::endl;
    }
    MTL::Buffer* sortedBuffer = device->newBuffer(nParticles * sizeof(bool), buffer->resourceOptions());
    bool success = false;
    MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
    cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);
    while (not success) {
        MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer(cmdDesc);
        MTL::ComputeCommandEncoder* computeCommandEncoder = commandBuffer->computeCommandEncoder();
        computeCommandEncoder->setBuffer(buffer, 0, 0);
        computeCommandEncoder->setBuffer(sortedBuffer, 0, 1);
        computeCommandEncoder->setBuffer(_cellArrayi, 0, 2);
        encodeCommand(computeCommandEncoder, _shuffleBool, nParticles);
        computeCommandEncoder->endEncoding();
        commandBuffer->commit();
        commandBuffer->waitUntilCompleted();
        success = commandBuffer->status() == MTL::CommandBufferStatusCompleted;
        if (!success) {
            std::cout << "Error, retrying" << std::endl;
        }
    }
    cmdDesc->release();
    buffer->release();
    return sortedBuffer;
}


void Compute::shuffleTreeLayer(MTL::Device* device, MTL::CommandQueue* commandQueue, int layer) {
    bool success = false;
    MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
    cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);

    while (not success) {
        MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer(cmdDesc);
        MTL::ComputeCommandEncoder* computeCommandEncoder = commandBuffer->computeCommandEncoder();
        computeCommandEncoder->setBuffer(_tree, 0, 0);
        computeCommandEncoder->setBuffer(_cellArrayj, 0, 1);
        computeCommandEncoder->setBuffer(_treeLevelBuffers[layer], 0, 2);
        encodeCommand(computeCommandEncoder, _shuffleTree, treeLevels[layer].size());
        computeCommandEncoder->endEncoding();
        commandBuffer->commit();
        commandBuffer->waitUntilCompleted();
        success = commandBuffer->status() == MTL::CommandBufferStatusCompleted;
        if (!success) {
            std::cout << "Error, retrying" << std::endl;
        }
    }
    cmdDesc->release();
}

void Compute::inverseCellParticles(MTL::Device *device, MTL::CommandQueue *commandQueue) {
    bool success = false;
    MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
    cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);

    while (not success) {
        MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer(cmdDesc);
        MTL::ComputeCommandEncoder* computeCommandEncoder = commandBuffer->computeCommandEncoder();
        computeCommandEncoder->setBuffer(_cellArrayi, 0, 0);
        computeCommandEncoder->setBuffer(_cellArrayj, 0, 1);
        encodeCommand(computeCommandEncoder, _inverseArray, nParticles);
        computeCommandEncoder->endEncoding();
        commandBuffer->commit();
        commandBuffer->waitUntilCompleted();
        success = commandBuffer->status() == MTL::CommandBufferStatusCompleted;
        if (!success) {
            std::cout << "Error, retrying" << std::endl;
        }
    }
    cmdDesc->release();
}

// -------------------------------- //
//                                  //
//          Octree Stuff            //
//                                  //
// -------------------------------- //

void Compute::updateOctreeBuffer(MTL::Device* device) {
    treeLevels.clear();
    treeLevels = treeLevelsTemp;
    
    long nodeValuesSize = nodeValues;
    if (nodeValuesSize > prevNodeValues) {
        // Multiply by a factor so we don't have to keep claiming large amounts of memory.
        nodeValuesSize *= EXTRA_MEMORY_MULTIPLIER;
        _multipoleExpansions->release();
        _localExpansion->release();
        _parentIndexes->release();
        _multipoleExpansions = device->newBuffer(nodeValuesSize * sizeof(Multipole), MTL::ResourceStorageModePrivate);
        _localExpansion = device->newBuffer(nodeValuesSize * sizeof(Local), MTL::ResourceStorageModePrivate);
        _parentIndexes = device->newBuffer(nodeValuesSize * sizeof(unsigned long), MTL::ResourceStorageModePrivate);
        prevNodeValues = nodeValuesSize;
    } else if (nodeValuesSize < prevNodeValues * RED_MEMEORY_MULTIPLIER) {
        nodeValuesSize *= EXTRA_MEMORY_MULTIPLIER;
        _multipoleExpansions->release();
        _localExpansion->release();
        _parentIndexes->release();
        _multipoleExpansions = device->newBuffer(nodeValuesSize * sizeof(Multipole), MTL::ResourceStorageModePrivate);
        _localExpansion = device->newBuffer(nodeValuesSize * sizeof(Local), MTL::ResourceStorageModePrivate);
        _parentIndexes = device->newBuffer(nodeValuesSize * sizeof(unsigned long), MTL::ResourceStorageModePrivate);
        prevNodeValues = nodeValuesSize;
    }

    _tree->release();
    _tree = _treeTmp;
    
    for (int i = 0; i < _treeLevelBuffers.size(); i++) {
        _treeLevelBuffers[i]->release();
    }
    _treeLevelBuffers.clear();
    long maxSizei = 0;
    long maxSizej = 0;
    for (int i = 0; i < treeLevels.size(); i++) {
        MTL::Buffer* treeLevelBuffer = device->newBuffer(treeLevels[i].size() * sizeof(int), MTL::ResourceStorageModeShared);
        if (i != treeLevels.size() - 1) {
            if (i % 2 == 0) {
                maxSizej = fmax(maxSizej, treeLevels[i].size());
            } else {
                maxSizei = fmax(maxSizei, treeLevels[i].size());
            }
        }
        writeDataToBuffer(treeLevelBuffer, treeLevels[i]);
        _treeLevelBuffers.push_back(treeLevelBuffer);
    }
    long gravDataSizei = maxSizei * sizeof(int) * MAX_UNCHECKED_POINTERS;
    long gravDataSizej = maxSizej * sizeof(int) * MAX_UNCHECKED_POINTERS;
    if (gravDataSizei > prevGravDataSizei) {
        gravDataSizei *= EXTRA_MEMORY_MULTIPLIER;
        _localGravi->release();
        _localGravi = device->newBuffer(gravDataSizei, MTL::ResourceStorageModePrivate);
        prevGravDataSizei = gravDataSizei;
    } else if (gravDataSizei < prevGravDataSizei * RED_MEMEORY_MULTIPLIER) {
        gravDataSizei *= EXTRA_MEMORY_MULTIPLIER;
        _localGravi->release();
        _localGravi = device->newBuffer(gravDataSizei, MTL::ResourceStorageModePrivate);
        prevGravDataSizei = gravDataSizei;
    }
    if (gravDataSizej > prevGravDataSizej) {
        gravDataSizej *= EXTRA_MEMORY_MULTIPLIER;
        _localGravj->release();
        _localGravj = device->newBuffer(gravDataSizej, MTL::ResourceStorageModePrivate);
        prevGravDataSizej = gravDataSizej;
    } else if (gravDataSizej < prevGravDataSizej * RED_MEMEORY_MULTIPLIER) {
        gravDataSizej *= EXTRA_MEMORY_MULTIPLIER;
        _localGravj->release();
        _localGravj = device->newBuffer(gravDataSizej, MTL::ResourceStorageModePrivate);
        prevGravDataSizej = gravDataSizej;
    }
    
    updating = true;
    std::thread([this](MTL::Device* _device) {
        this->updateOctreeData(_device);
    }, device).detach();
}

void Compute::saveState(MTL::Device* device, MTL::CommandQueue* commandQueue, std::string filename) {
    // Create shared buffers for all the private data we need to copy
    MTL::Buffer* velocities = device->newBuffer(sizeof(simd_float3) * nParticles, MTL::ResourceStorageModeShared);
    MTL::Buffer* densities = device->newBuffer(sizeof(float) * nParticles, MTL::ResourceStorageModeShared);
    MTL::Buffer* internalEnergy = device->newBuffer(sizeof(float) * nParticles, MTL::ResourceStorageModeShared);
    MTL::Buffer* masses = device->newBuffer(sizeof(float) * nParticles, MTL::ResourceStorageModeShared);
    MTL::Buffer* materialIds = device->newBuffer(sizeof(int) * nParticles, MTL::ResourceStorageModeShared);
    MTL::Buffer* pressures = device->newBuffer(sizeof(float) * nParticles, MTL::ResourceStorageModeShared);
    MTL::Buffer* smoothingLengths = device->newBuffer(sizeof(float) * nParticles, MTL::ResourceStorageModeShared);
    
    // Copy the data into them
    MTL::CommandBuffer* cmdBuffer = commandQueue->commandBuffer();
    MTL::BlitCommandEncoder* blitEncoder = cmdBuffer->blitCommandEncoder();
    blitEncoder->copyFromBuffer(_velocityBuffer, 0, velocities, 0, sizeof(simd_float3) * nParticles);
    blitEncoder->copyFromBuffer(densityBuffer, 0, densities, 0, sizeof(float) * nParticles);
    blitEncoder->copyFromBuffer(_internalEnergyBuffer, 0, internalEnergy, 0, sizeof(float) * nParticles);
    blitEncoder->copyFromBuffer(_massBuffer, 0, masses, 0, sizeof(float) * nParticles);
    blitEncoder->copyFromBuffer(materialIdBuffer, 0, materialIds, 0, sizeof(int) * nParticles);
    blitEncoder->copyFromBuffer(_pressureBuffer, 0, pressures, 0, sizeof(float) * nParticles);
    blitEncoder->copyFromBuffer(_smoothingLengthBuffer, 0, smoothingLengths, 0, sizeof(float) * nParticles);
    blitEncoder->endEncoding();
    cmdBuffer->commit();
    cmdBuffer->waitUntilCompleted();
    
    DataStruct data;
    simd_float3* poses = static_cast<simd_float3*>(positionBuffer->contents());
    data.positions.insert(data.positions.end(), &poses[0], &poses[nParticles]);
    simd_float3* vels = static_cast<simd_float3*>(velocities->contents());
    data.velocities.insert(data.velocities.end(), &vels[0], &vels[nParticles]);
    float* dens = static_cast<float*>(densities->contents());
    data.densities.insert(data.densities.end(), &dens[0], &dens[nParticles]);
    float* ie = static_cast<float*>(internalEnergy->contents());
    data.internalEnergy.insert(data.internalEnergy.end(), &ie[0], &ie[nParticles]);
    float* m = static_cast<float*>(masses->contents());
    data.masses.insert(data.masses.end(), &m[0], &m[nParticles]);
    int* mids = static_cast<int*>(materialIds->contents());
    data.materialIDs.insert(data.materialIDs.end(), &mids[0], &mids[nParticles]);
    float* p = static_cast<float*>(pressures->contents());
    data.pressures.insert(data.pressures.end(), &p[0], &p[nParticles]);
    float* h = static_cast<float*>(smoothingLengths->contents());
    data.smoothingLengths.insert(data.smoothingLengths.end(), &h[0], &h[nParticles]);
    
    writeHDFFile(filename, data);
    
    velocities->release();
    densities->release();
    internalEnergy->release();
    masses->release();
    materialIds->release();
    pressures->release();
    smoothingLengths->release();
}
