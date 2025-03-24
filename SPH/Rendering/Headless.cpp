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
    particles = new Particles(_device, compute->positionBuffer, compute->materialIdBuffer, compute->densityBuffer, compute->_massBuffer, compute->_smoothingLengthBuffer, compute->rhoGrads, compute->temperature, compute->_alpha, compute->nParticles);
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
    if (time > nextSnapshot) {
        int initialSnapshot = round(START_SNAPSHOT * (1000.f / SNAPSHOT_PERIOD_SECONDS));
        int currentFrame = round(nextSnapshot / SNAPSHOT_PERIOD_SECONDS);
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
        
        nextSnapshot += SNAPSHOT_PERIOD_SECONDS;
        return (nextSnapshot - SNAPSHOT_PERIOD_SECONDS) % 1000 == 0 and time > 500;
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
            compute->activatePass(pCmd1);
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
                    compute->densityPass(pCmd3, false);
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
            compute->accelerationStepPass(pCmd4);
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
            compute->activatePass(pCmd);
            compute->gravitationalPass(pCmd);
            if (_frame == 0) {
                // First frame loop density to get smooth gradients
                for (int j = 0; j < DENSITY_GRADIENT_SETTLING_ITTERATIONS; j++) {
                    compute->densityPass(pCmd, false);
                }
            }
            compute->densityPass(pCmd);
            compute->accelerationPass(pCmd);
            compute->accelerationStepPass(pCmd);
            compute->stepPass(pCmd);
            _frame ++;
        }
        pCmd->commit();
        pCmd->waitUntilCompleted();
    }
    
    bool save = takeSnapshot();
    compute->organisation(_device, _commandQueue);
    particles->updateBuffers(compute->positionBuffer, compute->materialIdBuffer, compute->densityBuffer, compute->_massBuffer, compute->_smoothingLengthBuffer, compute->rhoGrads, compute->temperature, compute->_alpha);
    if (save) {
        compute->saveState(_device, _commandQueue, std::string(SAVES_DIR) + "/state_" + std::to_string(START_SNAPSHOT + (nextSnapshot / 1000)) + ".hdf5");
    }
}

void Headless::run() {
    int i = 0;
    auto start = std::chrono::high_resolution_clock::now();
    auto const constStart = start;
    float simTimStart = getTime();
    const float constSimStart = simTimStart;
    int framesPerCout = HEADLESS_ITTERATION_RATE_REFRESH * fmax(1, STEPS_PER_FRAME);
    bool run = true;
    while (run) {
        step();
        if (i == HEADLESS_ITTERATION_RATE_REFRESH) {
            auto end = std::chrono::high_resolution_clock::now();
            float simTimeEnd = getTime();
            std::chrono::duration<double> elapsed = end - start;
            std::chrono::duration<double> totalElapsed = end - constStart;
            float simElapsed = simTimeEnd - simTimStart;
            float fps = (float)framesPerCout / elapsed.count();
            float secsPerSecs = simElapsed / elapsed.count();
            float avgSecsPerSecs = (simTimeEnd - constSimStart) / totalElapsed.count();
            float dtAvg = simElapsed / framesPerCout;
            std::cout << "Steps Per Second: " << std::to_string(fps);
            std::cout << "\nSimulations seconds per second: " << std::to_string(secsPerSecs);
            std::cout << "\nAverage simulations seconds per second: " << std::to_string(avgSecsPerSecs);
            std::cout << "\nAverage dt: " << std::to_string(dtAvg);
            std::cout << "\nSimulation time: " << std::to_string(START_SNAPSHOT * 1000 + simTimeEnd) << "\n" << std::endl;
            i = 0;
            start = std::chrono::high_resolution_clock::now();
            simTimStart = simTimeEnd;
        }
        i++;
    }
}
