//
//  Snapshotter.cpp
//  SPH
//
//  Created by Charlie Close on 09/03/2025.
//

#include "Snapshotter.hpp"
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

Snapshotter::Snapshotter(MTL::Device* device, Particles* particles, simd_float3 position, float pitch, float yaw, std::string snapshotFolder) {
    _device = device;
    _commandQueue = _device->newCommandQueue();
    _particles = particles;
    _camera = new Camera(position, pitch, yaw);
    _snapshotFolder = snapshotFolder;
    
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
    *cameraData = _camera->getMatrix();
    _cameraPosBuffer = _device->newBuffer( sizeof(simd::float3), MTL::ResourceStorageModeShared );
    simd::float3* cameraPos = static_cast<simd::float3*>(_cameraPosBuffer->contents());
    *cameraPos = _camera->getPosition();
    
    createRenderPassDescriptor();
    
    if (!std::filesystem::is_directory(SNAPSHOT_DIR + std::string("/") + snapshotFolder)) {
        std::filesystem::create_directory(SNAPSHOT_DIR + std::string("/") + snapshotFolder);
    }
}

Snapshotter::~Snapshotter()
{
    _commandQueue->release();
}

void Snapshotter::createRenderPassDescriptor() {
    renderPassDescriptor = MTL::RenderPassDescriptor::renderPassDescriptor();
    MTL::RenderPassDepthAttachmentDescriptor* depthAttachementDescriptor = renderPassDescriptor->depthAttachment();
    MTL::RenderPassColorAttachmentDescriptorArray* colorAttachmentDescriptorArray = renderPassDescriptor->colorAttachments();
    MTL::RenderPassColorAttachmentDescriptor* colorAttachmentDescriptor = colorAttachmentDescriptorArray->object(0);
    depthAttachementDescriptor->setTexture(depthTexture);
    depthAttachementDescriptor->setLoadAction(MTL::LoadActionClear);
    colorAttachmentDescriptor->setTexture(colourTexture);
    colorAttachmentDescriptor->setLoadAction(MTL::LoadActionClear);
}

void Snapshotter::takeSnapshot(int n) {
    while (loading) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    loading = true;
    while (true) {
        MTL::CommandBufferDescriptor* cmdDesc = MTL::CommandBufferDescriptor::alloc();
        cmdDesc->setErrorOptions(MTL::CommandBufferErrorOptionEncoderExecutionStatus);
        MTL::CommandBuffer* pCmd = _commandQueue->commandBuffer(cmdDesc);
        _particles->zeroVis(pCmd);
        _particles->calculateNormals(pCmd, _cameraPosBuffer);
        _particles->setIndrBuff(pCmd);
        MTL::RenderCommandEncoder* cmdEncoder = pCmd->renderCommandEncoder(renderPassDescriptor);
        _particles->draw(cmdEncoder, _cameraDataBuffer, _cameraPosBuffer);
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
        lodepng::save_file(encodedData, std::string(SNAPSHOT_DIR)  + std::string("/") + _snapshotFolder + "/snapshot_" + std::to_string(n) + ".png");
    }, n).detach();
}
