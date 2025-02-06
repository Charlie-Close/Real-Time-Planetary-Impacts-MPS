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

// -------------------------------- //
//                                  //
//          Constructor             //
//                                  //
// -------------------------------- //

Compute::Compute(MTL::Device* device) {
//    DataStruct data = readHDFFile("earth_impact_0000.hdf5");
//    DataStruct data = readHDFFile("demo_impact_n60.hdf5");
    DataStruct data = readHDFFile("demo_impact_n50.hdf5");
//    DataStruct data = readHDFFile("demo_impact_n40.hdf5");
    
    MTL::CommandQueue* commandQueue = device->newCommandQueue();
        
    nParticles = (int)data.positions.size();
    cellSize = .2;
    cellsPerDim = pow(2, 7);
    
    std::cout << nParticles << std::endl;

    buildShaders(device);
    buildBuffers(device, commandQueue);
    loadInitialConditions(device, commandQueue, data);
}

// -------------------------------- //
//                                  //
//   Initial Condition Generation   //z
//                                  //
// -------------------------------- //

void Compute::loadInitialConditions(MTL::Device* device, MTL::CommandQueue* commandQueue, DataStruct data) {
    ANEOSTable forTable = loadANEOSDataFromFile("ANEOS_forsterite_S19.txt", 1024);
    ANEOSTable feTable = loadANEOSDataFromFile("ANEOS_Fe85Si15_S20.txt", 1024);
    
    _forsterite = createRG32FloatTexture(device, commandQueue, forTable);
    _Fe85Si15 = createRG32FloatTexture(device, commandQueue, feTable);
    
    // Write all the data to the buffers
    writeDataToBuffer(positionBuffer, data.positions);
    writeDataToPrivateBuffer(device, commandQueue, _velocityBuffer, data.velocities);
    writeDataToPrivateBuffer(device, commandQueue, _densityBuffer, data.densities);
    writeDataToPrivateBuffer(device, commandQueue, _internalEnergyBuffer, data.internalEnergy);
    writeDataToPrivateBuffer(device, commandQueue, _massBuffer, data.masses);
    writeDataToPrivateBuffer(device, commandQueue, _pressureBuffer, data.pressures);
    writeDataToPrivateBuffer(device, commandQueue, _smoothingLengthBuffer, data.smoothingLengths);
    writeDataToPrivateBuffer(device, commandQueue, _materialIdBuffer, data.materialIDs);
    writeDataToPrivateBuffer(device, commandQueue, _nParticles, nParticles);
    for (uint i = 0; i < 4; i++) {
        writeDataToPrivateBuffer(device, commandQueue, _ittr[i], i);
    }
    writeDataToPrivateBuffer(device, commandQueue, _nBlocks, nBlocks);
    
    bool isAlive[nParticles];
    for (int i = 0; i < nParticles; i++) {
        isAlive[i] = true;
    }
    writeDataToPrivateBuffer(device, commandQueue, _isAlive, isAlive, nParticles);
    writeDataToBuffer(_debug, 0);
    writeDataToBuffer(_false, false);
    writeDataToBuffer(_true, true);

    
    writeDataToPrivateBuffer(device, commandQueue, _cellSize, cellSize);
    writeDataToPrivateBuffer(device, commandQueue, _cellsPerDim, cellsPerDim);
    
    updateOctreeData(device);
    updateOctreeBuffer(device);
}

void Compute::updateOctreeBuffer(MTL::Device* device) {
    if (updating) {
        return;
    }
    
    treeLevels = treeLevelsTemp;
    
    long treeDataSize = nodeValues;
    if (treeDataSize > prevTreeValuesSize) {
        treeDataSize *= 1.2f;
        _multipoleExpansions->release();
        _localExpansion->release();
        _parentIndexes->release();
        _multipoleExpansions = device->newBuffer(treeDataSize * sizeof(Multipole), MTL::ResourceStorageModePrivate);
        _localExpansion = device->newBuffer(treeDataSize * sizeof(Local), MTL::ResourceStorageModePrivate);
        _parentIndexes = device->newBuffer(treeDataSize * sizeof(uint), MTL::ResourceStorageModePrivate);
        prevTreeValuesSize = treeDataSize;
    }

    _tree->release();
    _tree = _treeTmp;
    
    _treeLevelBuffers.clear();
    long maxSize = 0;
    for (int i = 0; i < treeLevels.size(); i++) {
        MTL::Buffer* treeLevelBuffer = device->newBuffer(treeLevels[i].size() * sizeof(int), MTL::ResourceStorageModeShared);
        if (i != treeLevels.size() - 1) {
            maxSize = fmax(maxSize, treeLevels[i].size());
        }
        writeDataToBuffer(treeLevelBuffer, treeLevels[i]);
        _treeLevelBuffers.push_back(treeLevelBuffer);
    }
    
    long gravDataSize = maxSize * sizeof(int) * 8192;
    if (gravDataSize > prevGravDataSize) {
        gravDataSize *= 1.2;
        _localGravi->release();
        _localGravj->release();
        _localGravi = device->newBuffer(gravDataSize, MTL::ResourceStorageModePrivate);
        _localGravj = device->newBuffer(gravDataSize, MTL::ResourceStorageModePrivate);
        prevGravDataSize = gravDataSize;
    }
    
    updating = true;
    std::thread([this](MTL::Device* _device) {
        this->updateOctreeData(_device);
    }, device).detach();
}

void Compute::updateOctreeBufferTest(MTL::Device* device) {
//    if (updating) {
//        return;
//    }
//    
//    treeLevels = treeLevelsTemp;
//    
//    long treeDataSize = nodeValues;
//    if (treeDataSize > prevTreeValuesSize) {
//        treeDataSize *= 1.5f;
//        _multipoleExpansions->release();
//        _localExpansion->release();
//        _parentIndexes->release();
//        _multipoleExpansions = device->newBuffer(treeDataSize * sizeof(Multipole), MTL::ResourceStorageModePrivate);
//        _localExpansion = device->newBuffer(treeDataSize * sizeof(Local), MTL::ResourceStorageModePrivate);
//        _parentIndexes = device->newBuffer(treeDataSize * sizeof(uint), MTL::ResourceStorageModePrivate);
//        prevTreeValuesSize = treeDataSize;
//    }
//
//    _tree->release();
//    _tree = _treeTmp;
//    
//    _treeLevelBuffers.clear();
//    long maxSize = 0;
//    for (int i = 0; i < treeLevels.size(); i++) {
//        MTL::Buffer* treeLevelBuffer = device->newBuffer(treeLevels[i].size() * sizeof(int), MTL::ResourceStorageModeShared);
//        if (i != treeLevels.size() - 1) {
//            maxSize = fmax(maxSize, treeLevels[i].size());
//        }
//        writeDataToBuffer(treeLevelBuffer, treeLevels[i]);
//        _treeLevelBuffers.push_back(treeLevelBuffer);
//    }
//    
//    long gravDataSize = maxSize * sizeof(int) * 8192;
//    if (gravDataSize > prevGravDataSize) {
//        gravDataSize *= 1.5;
//        _localGravi->release();
//        _localGravj->release();
//        _localGravi = device->newBuffer(gravDataSize, MTL::ResourceStorageModePrivate);
//        _localGravj = device->newBuffer(gravDataSize, MTL::ResourceStorageModePrivate);
//        prevGravDataSize = gravDataSize;
//    }
//    
//    updating = true;
//    std::thread([this](MTL::Device* _device) {
//        this->updateOctreeData(_device);
//    }, device).detach();
}


void Compute::updateOctreeData(MTL::Device* device) {
    simd_float3* positions = static_cast<simd_float3*>(positionBuffer->contents());
    bool* isAlive = static_cast<bool*>(_isAlive->contents());
    buildOctree(positions, isAlive, nParticles, octreeData, treeLevelsTemp, nodeValues, 8);
    _treeTmp = device->newBuffer(octreeData.size() * sizeof(int), MTL::ResourceStorageModeShared);
    writeDataToBuffer(_treeTmp, octreeData);
    updating = false;
}

// -------------------------------- //
//                                  //
//      Metal Object Builders       //
//                                  //
// -------------------------------- //

void Compute::buildBuffers(MTL::Device* device, MTL::CommandQueue *commandQueue) {
    nBlocks = (nParticles / 256) + 1;
    std::cout << "N BLOCK " << nBlocks << std::endl;
    positionBuffer = device->newBuffer(nParticles * sizeof(simd_int3), MTL::ResourceStorageModeShared);
    _velocityBuffer = device->newBuffer(nParticles * sizeof(simd_int3), MTL::ResourceStorageModePrivate);
    _accelerationBuffer = device->newBuffer(nParticles * sizeof(simd_int3), MTL::ResourceStorageModePrivate);
    _densityBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _internalEnergyBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _massBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _pressureBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _smoothingLengthBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _materialIdBuffer = device->newBuffer(nParticles * sizeof(int), MTL::ResourceStorageModePrivate);
    _gradientTermsBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _speedOfSoundBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _dInternalEnergyBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _potentialGradientBuffer = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _balsara = device->newBuffer(nParticles * sizeof(float), MTL::ResourceStorageModePrivate);
    _cellArrayi = device->newBuffer(nParticles * sizeof(simd_uint2), MTL::ResourceStorageModePrivate);
    _cellArrayj = device->newBuffer(nParticles * sizeof(simd_uint2), MTL::ResourceStorageModePrivate);
    _bucketHist = device->newBuffer(nBlocks * 256 * sizeof(uint), MTL::ResourceStorageModePrivate);
    _bucketOffset = device->newBuffer(256 * sizeof(uint), MTL::ResourceStorageModePrivate);
    _particleOffset = device->newBuffer(nParticles * sizeof(uint), MTL::ResourceStorageModePrivate);
    for (int i = 0; i < 4; i++) {
        _ittr[i] = device->newBuffer(sizeof(uint), MTL::ResourceStorageModePrivate);
    }
    _nParticles = device->newBuffer(sizeof(int), MTL::ResourceStorageModePrivate);
    _nBlocks = device->newBuffer(sizeof(int), MTL::ResourceStorageModePrivate);
    _cellStart = device->newBuffer(cellsPerDim * cellsPerDim * cellsPerDim * sizeof(int), MTL::ResourceStorageModePrivate);
    _cellEnd = device->newBuffer(cellsPerDim * cellsPerDim * cellsPerDim * sizeof(int), MTL::ResourceStorageModePrivate);
    _cellSize = device->newBuffer(sizeof(float), MTL::ResourceStorageModePrivate);
    _cellsPerDim = device->newBuffer(sizeof(int), MTL::ResourceStorageModePrivate);
    _leafPointers = device->newBuffer(nParticles * sizeof(int), MTL::ResourceStorageModePrivate);
    _isAlive = device->newBuffer(sizeof(bool) * nParticles, MTL::ResourceStorageModeShared);
    _debug = device->newBuffer(sizeof(int), MTL::ResourceStorageModeShared);
    _false = device->newBuffer(sizeof(bool), MTL::ResourceStorageModeShared);
    _true = device->newBuffer(sizeof(bool), MTL::ResourceStorageModeShared);
    _dt = device->newBuffer(sizeof(float), MTL::ResourceStorageModePrivate);
    _gravAbs = device->newBuffer(sizeof(float) * nParticles, MTL::ResourceStorageModePrivate);
    _tree = nullptr;

}

void Compute::buildShaders(MTL::Device* device) {
    NS::Error** error = nil;

    // Load the shader files with a .metal file extension in the project
    MTL::Library* defaultLibrary = device->newDefaultLibrary();
        
    MTL::Function* densityFunction = defaultLibrary->newFunction(NS::String::string("density", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* accelerationFunction = defaultLibrary->newFunction(NS::String::string("acceleration", NS::StringEncoding::UTF8StringEncoding));
//    MTL::Function* gravityFunction = defaultLibrary->newFunction(NS::String::string("gravity", NS::StringEncoding::UTF8StringEncoding));
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


    
    // Create a compute pipeline state object.
    _densityPSO = device->newComputePipelineState(densityFunction, error);
    _accelerationPSO = device->newComputePipelineState(accelerationFunction, error);
//    _gravityPSO = device->newComputePipelineState(gravityFunction, error);
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
    computeEncoder->setBuffer(_densityBuffer, 0, 2);
    computeEncoder->setBuffer(_internalEnergyBuffer, 0, 3);
    computeEncoder->setBuffer(_massBuffer, 0, 4);
    computeEncoder->setBuffer(_pressureBuffer, 0, 5);
    computeEncoder->setBuffer(_smoothingLengthBuffer, 0, 6);
    computeEncoder->setBuffer(_materialIdBuffer, 0, 7);
    computeEncoder->setBuffer(_gradientTermsBuffer, 0, 8);
    computeEncoder->setBuffer(_speedOfSoundBuffer, 0, 9);
    computeEncoder->setBuffer(_potentialGradientBuffer, 0, 10);
    computeEncoder->setBuffer(_balsara, 0, 11);
    computeEncoder->setBuffer(_cellArrayi, 0, 12);
    computeEncoder->setBuffer(_cellStart, 0, 13);
    computeEncoder->setBuffer(_cellEnd, 0, 14);
    computeEncoder->setBuffer(_cellSize, 0, 15);
    computeEncoder->setBuffer(_cellsPerDim, 0, 16);
    computeEncoder->setBuffer(_isAlive, 0, 17);
    computeEncoder->setBuffer(_dt, 0, 18);
    computeEncoder->setTexture(_forsterite, 0);
    computeEncoder->setTexture(_Fe85Si15, 1);
    
    encodeCommand(computeEncoder, _densityPSO, nParticles);

    //  End the compute pass.
    computeEncoder->endEncoding();
}

void Compute::accelerationPass(MTL::CommandBuffer* commandBuffer) {
    MTL::ComputeCommandEncoder* computeEncoder = commandBuffer->computeCommandEncoder();
    assert(computeEncoder != nil);
    
    computeEncoder->setBuffer(positionBuffer, 0, 0);
    computeEncoder->setBuffer(_velocityBuffer, 0, 1);
    computeEncoder->setBuffer(_accelerationBuffer, 0, 2);
    computeEncoder->setBuffer(_densityBuffer, 0, 3);
    computeEncoder->setBuffer(_internalEnergyBuffer, 0, 4);
    computeEncoder->setBuffer(_massBuffer, 0, 5);
    computeEncoder->setBuffer(_pressureBuffer, 0, 6);
    computeEncoder->setBuffer(_smoothingLengthBuffer, 0, 7);
    computeEncoder->setBuffer(_gradientTermsBuffer, 0, 8);
    computeEncoder->setBuffer(_speedOfSoundBuffer, 0, 9);
    computeEncoder->setBuffer(_dInternalEnergyBuffer, 0, 10);
    computeEncoder->setBuffer(_potentialGradientBuffer, 0, 11);
    computeEncoder->setBuffer(_balsara, 0, 12);
    computeEncoder->setBuffer(_cellArrayi, 0, 13);
    computeEncoder->setBuffer(_cellStart, 0, 14);
    computeEncoder->setBuffer(_cellEnd, 0, 15);
    computeEncoder->setBuffer(_cellSize, 0, 16);
    computeEncoder->setBuffer(_cellsPerDim, 0, 17);
    computeEncoder->setBuffer(_isAlive, 0, 18);
    computeEncoder->setBuffer(_dt, 0, 19);
    encodeCommand(computeEncoder, _accelerationPSO, nParticles);

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
    computeEncoder->setBuffer(_internalEnergyBuffer, 0, 3);
    computeEncoder->setBuffer(_dInternalEnergyBuffer, 0, 4);
    computeEncoder->setBuffer(_isAlive, 0, 5);
    computeEncoder->setBuffer(_dt, 0, 6);
    encodeCommand(computeEncoder, _mStepPSO, nParticles);

    //  End the compute pass.
    computeEncoder->endEncoding();
}

void Compute::sort(MTL::CommandBuffer *commandBuffer) {
    // Calcualte particle indexes and write to cell data
    
    MTL::ComputeCommandEncoder* computeEncoder = commandBuffer->computeCommandEncoder();
    computeEncoder->setBuffer(positionBuffer, 0, 0);
    computeEncoder->setBuffer(_cellSize, 0, 1);
    computeEncoder->setBuffer(_cellsPerDim, 0, 2);
    computeEncoder->setBuffer(_cellArrayi, 0, 3);
    encodeCommand(computeEncoder, _hashPSO, nParticles);
    computeEncoder->endEncoding();
    
    for (int i = 0; i < 2; i++) {
        
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
        encodeCommand(cEncScan, _scanPSO, 256);
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
    for (int i = _treeLevelBuffers.size() - 1; i >= 0; i--) {
        MTL::ComputeCommandEncoder* computeEncoder = commandBuffer->computeCommandEncoder();
        assert(computeEncoder != nil);
        
        computeEncoder->setBuffer(positionBuffer, 0, 0);
        computeEncoder->setBuffer(_massBuffer, 0, 1);
        computeEncoder->setBuffer(_tree, 0, 2);
        computeEncoder->setBuffer(_multipoleExpansions, 0, 3);
        computeEncoder->setBuffer(_localExpansion, 0, 4);
        computeEncoder->setBuffer(_treeLevelBuffers[i], 0, 5);
        computeEncoder->setBuffer(_parentIndexes, 0, 6);
        computeEncoder->setBuffer(_gravAbs, 0, 7);

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
        computeEncoder->setBuffer(_tree, 0, 3);
        computeEncoder->setBuffer(_multipoleExpansions, 0, 4);
        computeEncoder->setBuffer(_localExpansion, 0, 5);
        computeEncoder->setBuffer(_treeLevelBuffers[i], 0, 6);
        computeEncoder->setBuffer(_parentIndexes, 0, 7);
        computeEncoder->setBuffer(i % 2 == 0 ? _localGravi : _localGravj, 0, 8);
        computeEncoder->setBuffer(i % 2 == 0 ? _localGravj : _localGravi, 0, 9);
        computeEncoder->setBuffer(_gravAbs, 0, 10);


        encodeCommand(computeEncoder, _downTreePSO, treeLevels[i].size());

        //  End the compute pass.
        computeEncoder->endEncoding();
    }
}

void Compute::encodeCommand(MTL::ComputeCommandEncoder* computeEncoder, MTL::ComputePipelineState* command, int size) {
    MTL::Size gridSize = MTL::Size(size, 1, 1);

    computeEncoder->setComputePipelineState(command);
    // Calculate a threadgroup size.
    NS::UInteger threadGroupSize = command->maxTotalThreadsPerThreadgroup();
    if (threadGroupSize > size)
    {
        threadGroupSize = size;
    }
    MTL::Size threadgroupSize = MTL::Size(threadGroupSize, 1, 1);

    // Encode the compute command.
    computeEncoder->dispatchThreads(gridSize, threadgroupSize);
}
