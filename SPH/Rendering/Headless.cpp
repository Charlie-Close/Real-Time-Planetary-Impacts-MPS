//
//  Headless.cpp
//  SPH
//
//  Created by Charlie Close on 10/02/2025.
//

#include "Headless.hpp"
#include "Buffers.hpp"
#include "lodepng.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <iostream>
#include <filesystem>

// -------------------------------- //
//                                  //
//          Constructor             //
//                                  //
// -------------------------------- //

Headless::Headless() {
    _device = MTL::CreateSystemDefaultDevice();
    _commandQueue = _device->newCommandQueue();
    compute = new Compute(_device);
    particles = new Particles(_device, compute->positionBuffer, compute->materialIdBuffer, compute->densityBuffer, compute->_massBuffer, compute->_smoothingLengthBuffer, compute->rhoGrads, compute->temperature, compute->nParticles);
    camera = new Camera();
    
    _cameraDataBuffer = _device->newBuffer( sizeof(simd::float4x4), MTL::ResourceStorageModeShared );
    simd::float4x4* cameraData = static_cast<simd::float4x4*>(_cameraDataBuffer->contents());
    *cameraData = camera->getMatrix();
    _cameraPosBuffer = _device->newBuffer( sizeof(simd::float3), MTL::ResourceStorageModeShared );
    simd::float3* cameraPos = static_cast<simd::float3*>(_cameraPosBuffer->contents());
    *cameraPos = camera->getPosition();
        
    if (!std::filesystem::is_directory(SAVES_DIR)) {
        std::filesystem::create_directory(SAVES_DIR);
    }
    if (!std::filesystem::is_directory(SNAPSHOT_DIR)) {
        std::filesystem::create_directory(SNAPSHOT_DIR);
    }
//    MTL::Device* device, Particles* particles, simd_float3 position, float pitch, float yaw, std::string snapshotFolder
    snapshotters = new Snapshotter*[N_SNAPSHOTTERS];
    simd::float3 positions[] = POSITIONS;
    float pithces[] = PITCHES;
    float yaws[] = YAWS;
    for (int i = 0; i < N_SNAPSHOTTERS; i++) {
        snapshotters[i] = new Snapshotter(_device, particles, positions[i], pithces[i], yaws[i], "snapshotter_" + std::to_string(i));
        
    }
}

Headless::~Headless()
{
    _commandQueue->release();
    _device->release();
    delete compute;
}

float Headless::getTime() {
    return compute->getTime();
}

bool Headless::takeSnapshot() {
    float time = getTime();
    const int step = 50;
    if (time > nextSnapshot) {
        int initialSnapshot = round(START_SNAPSHOT * (1000.f / step));
        int currentFrame = round(nextSnapshot / step);
        while (true) {
            MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
            cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);
            MTL::CommandBuffer* pCmd = _commandQueue->commandBuffer(cmdDesc);
            particles->zeroVis(pCmd);
            particles->calculateNormals(pCmd, _cameraPosBuffer);
            particles->setIndrBuff(pCmd);
            particles->drawShadowMap(pCmd);
            pCmd->commit();
            pCmd->waitUntilCompleted();
            if (pCmd->status() != MTL::CommandBufferStatusCompleted) {
                std::cout << "Snapshot error, retrying..." << std::endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            } else {
                break;
            }
        }
        for (int i = 0; i < N_SNAPSHOTTERS; i++) {
            snapshotters[i]->takeSnapshot(initialSnapshot + currentFrame);
        }
        
        nextSnapshot += step;
        return (nextSnapshot - step) % 1000 == 0;
    }
    return false;
}

void Headless::step()
{
    
    if (STEPS_PER_FRAME == 0) {
        MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
        cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);
        while (true) {
            MTL::CommandBuffer* pCmd1 = _commandQueue->commandBuffer(cmdDesc);
            compute->sort(pCmd1);
            pCmd1->commit();
            pCmd1->waitUntilCompleted();
            if (pCmd1->status() != MTL::CommandBufferStatusCompleted) {
                std::cout << "Sorting error, retrying..." << std::endl;
            } else {
                break;
            }
        }
        while (true) {
            MTL::CommandBuffer* pCmd2 = _commandQueue->commandBuffer(cmdDesc);
            compute->gravitationalPass(pCmd2);
            pCmd2->commit();
            pCmd2->waitUntilCompleted();
            if (pCmd2->status() != MTL::CommandBufferStatusCompleted) {
                std::cout << "Gravity error, retrying..." << std::endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            } else {
                break;
            }
        }
        if (_frame == 0) {
            for (int i = 0; i < DENSITY_GRADIENT_SETTLING_ITTERATIONS; i++) {
                while (true) {
                    MTL::CommandBuffer* pCmd3 = _commandQueue->commandBuffer(cmdDesc);
                    compute->densityPass(pCmd3);
                    pCmd3->commit();
                    pCmd3->waitUntilCompleted();
                    if (pCmd3->status() != MTL::CommandBufferStatusCompleted) {
                        std::cout << "Density error, retrying..." << std::endl;
                        std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    } else {
                        break;
                    }
                }
            }
        }
        while (true) {
            MTL::CommandBuffer* pCmd3 = _commandQueue->commandBuffer(cmdDesc);
            compute->densityPass(pCmd3);
            pCmd3->commit();
            pCmd3->waitUntilCompleted();
            if (pCmd3->status() != MTL::CommandBufferStatusCompleted) {
                std::cout << "Density error, retrying..." << std::endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            } else {
                break;
            }
        }
        while (true) {
            MTL::CommandBuffer* pCmd4 = _commandQueue->commandBuffer(cmdDesc);
            compute->accelerationPass(pCmd4);
            compute->stepPass(pCmd4);
            pCmd4->commit();
            pCmd4->waitUntilCompleted();
            if (pCmd4->status() != MTL::CommandBufferStatusCompleted) {
                std::cout << "Acceleration error, aborting..." << std::endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            } else {
                break;
            }
        }
        cmdDesc->release();
        _frame ++;
    } else {
        
        MTL::CommandBuffer* pCmd = _commandQueue->commandBuffer();
        
        for (int i = 0; i < STEPS_PER_FRAME; i++) {
            compute->sort(pCmd);
            compute->gravitationalPass(pCmd);
            if (_frame == 0) {
                // First frame loop density to get smooth gradients
                for (int j = 0; j < DENSITY_GRADIENT_SETTLING_ITTERATIONS; j++) {
                    compute->densityPass(pCmd);
                }
            }
            compute->densityPass(pCmd);
            compute->accelerationPass(pCmd);
            compute->stepPass(pCmd);
            _frame ++;
        }
        pCmd->commit();
        pCmd->waitUntilCompleted();
    }
    
    bool save = takeSnapshot();
    compute->organisation(_device, _commandQueue);
    particles->updateBuffers(compute->positionBuffer, compute->materialIdBuffer, compute->densityBuffer, compute->_massBuffer, compute->_smoothingLengthBuffer, compute->rhoGrads, compute->temperature);
    if (save) {
        compute->saveState(_device, _commandQueue, std::string(SAVES_DIR) + "/state_" + std::to_string(START_SNAPSHOT + (nextSnapshot / 1000)) + ".hdf5");
    }
}


