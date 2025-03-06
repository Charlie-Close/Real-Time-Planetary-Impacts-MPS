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
    particles = new Particles(_device, compute->positionBuffer, compute->materialIdBuffer, compute->densityBuffer, compute->_massBuffer, compute->_smoothingLengthBuffer, compute->rhoGrads, compute->_cellStart, compute->_cellEnd, compute->_cellArrayi, compute->nParticles);
    camera = new Camera();
    
    // Color descriptor
    MTL::TextureDescriptor* colorDesc = MTL::TextureDescriptor::texture2DDescriptor(
        MTL::PixelFormatRGBA8Unorm,  // typical color format
        SNAPSHOT_RESOLUTION,
        SNAPSHOT_RESOLUTION,
        false
    );
    colorDesc->setUsage(MTL::TextureUsageRenderTarget | MTL::TextureUsageShaderRead);
    colourTexture = _device->newTexture(colorDesc);

    // Depth descriptor
    MTL::TextureDescriptor* depthDesc = MTL::TextureDescriptor::texture2DDescriptor(
        MTL::PixelFormatDepth32Float,  // typical depth format
        SNAPSHOT_RESOLUTION,
        SNAPSHOT_RESOLUTION,
        false
    );
    depthDesc->setUsage(MTL::TextureUsageRenderTarget);
    depthTexture = _device->newTexture(depthDesc);
    
    _cameraDataBuffer = _device->newBuffer( sizeof(simd::float4x4), MTL::ResourceStorageModeShared );
    simd::float4x4* cameraData = static_cast<simd::float4x4*>(_cameraDataBuffer->contents());
    *cameraData = camera->getMatrix();
    _cameraPosBuffer = _device->newBuffer( sizeof(simd::float3), MTL::ResourceStorageModeShared );
    simd::float3* cameraPos = static_cast<simd::float3*>(_cameraPosBuffer->contents());
    *cameraPos = camera->getPosition();
    
    createRenderPassDescriptor();
    
    if (!std::filesystem::is_directory(SAVES_DIR)) {
        std::filesystem::create_directory(SAVES_DIR);
    }
    if (!std::filesystem::is_directory(SNAPSHOT_DIR)) {
        std::filesystem::create_directory(SNAPSHOT_DIR);
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

std::pair<MTL::Texture*, MTL::Texture*> Headless::createTextures() {
    // Color descriptor
    MTL::TextureDescriptor* colorDesc = MTL::TextureDescriptor::texture2DDescriptor(
        MTL::PixelFormatRGBA8Unorm,  // typical color format
        SNAPSHOT_RESOLUTION,
        SNAPSHOT_RESOLUTION,
        false
    );
    colorDesc->setUsage(MTL::TextureUsageRenderTarget | MTL::TextureUsageShaderRead);
    colourTexture = _device->newTexture(colorDesc);

    // Depth descriptor
    MTL::TextureDescriptor* depthDesc = MTL::TextureDescriptor::texture2DDescriptor(
        MTL::PixelFormatDepth32Float,  // typical depth format
        SNAPSHOT_RESOLUTION,
        SNAPSHOT_RESOLUTION,
        false
    );
    depthDesc->setUsage(MTL::TextureUsageRenderTarget);
    depthTexture = _device->newTexture(depthDesc);
    
    return { colourTexture, depthTexture };
}

void Headless::createRenderPassDescriptor() {
    renderPassDescriptor = MTL::RenderPassDescriptor::renderPassDescriptor();
    MTL::RenderPassDepthAttachmentDescriptor* depthAttachementDescriptor = renderPassDescriptor->depthAttachment();
    MTL::RenderPassColorAttachmentDescriptorArray* colorAttachmentDescriptorArray = renderPassDescriptor->colorAttachments();
    MTL::RenderPassColorAttachmentDescriptor* colorAttachmentDescriptor = colorAttachmentDescriptorArray->object(0);
    depthAttachementDescriptor->setTexture(depthTexture);
    depthAttachementDescriptor->setLoadAction(MTL::LoadActionClear);
    colorAttachmentDescriptor->setTexture(colourTexture);
    colorAttachmentDescriptor->setLoadAction(MTL::LoadActionClear);
}

bool Headless::takeSnapshot() {
    float time = getTime();
    const int step = 50;
    if (time > nextSnapshot) {
        int initialSnapshot = round(START_SNAPSHOT * (1000.f / step));
        int currentFrame = round(nextSnapshot / step);
        while (loading) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        loading = true;
        while (true) {
            MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
            cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);
            MTL::CommandBuffer* pCmd = _commandQueue->commandBuffer(cmdDesc);
            particles->zeroVis(pCmd);
            particles->calculateNormals(pCmd, _cameraPosBuffer);
            particles->calculateNormals(pCmd, _cameraPosBuffer);
            particles->setIndrBuff(pCmd);
            particles->drawShadowMap(pCmd);
            MTL::RenderCommandEncoder* cmdEncoder = pCmd->renderCommandEncoder(renderPassDescriptor);
            particles->draw(cmdEncoder, _cameraDataBuffer, _cameraPosBuffer);
            cmdEncoder->endEncoding();
            pCmd->commit();
            pCmd->waitUntilCompleted();
            if (pCmd->status() != MTL::CommandBufferStatusCompleted) {
                std::cout << "Snapshot error, retrying..." << std::endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            } else {
                break;
            }
        }
        
        std::thread([this](int n) {
            int pixelByteCount = 4 * sizeof(std::uint8_t);
            int imageBytesPerRow = SNAPSHOT_RESOLUTION * pixelByteCount;
            std::vector<std::uint8_t> dataIn(4 * SNAPSHOT_RESOLUTION * SNAPSHOT_RESOLUTION);
            colourTexture->getBytes(dataIn.data(), imageBytesPerRow, MTL::Region::Make2D(0, 0, SNAPSHOT_RESOLUTION, SNAPSHOT_RESOLUTION), 0);
            this->loading = false;
            std::vector<std::uint8_t> encodedData;
            lodepng::encode(encodedData, dataIn, SNAPSHOT_RESOLUTION, SNAPSHOT_RESOLUTION);
            lodepng::save_file(encodedData, std::string(SNAPSHOT_DIR) + "/snapshot_" + std::to_string(n) + ".png");
        }, initialSnapshot + currentFrame).detach();
        
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
    particles->updateBuffers(compute->positionBuffer, compute->materialIdBuffer, compute->densityBuffer, compute->_massBuffer, compute->_smoothingLengthBuffer, compute->rhoGrads);
    if (save) {
        compute->saveState(_device, _commandQueue, std::string(SAVES_DIR) + "/state_" + std::to_string(START_SNAPSHOT + (nextSnapshot / 1000)) + ".hdf5");
    }
}


