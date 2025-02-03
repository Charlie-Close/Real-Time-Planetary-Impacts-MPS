//
//  Buffers.cpp
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include "Buffers.hpp"

void copyDataToBuffer(MTL::CommandQueue* commandQueue, MTL::Buffer *privateBuffer, MTL::Buffer *sharedBuffer, unsigned long bufferLength) {
    MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer();
    
    MTL::BlitCommandEncoder* blitCommandEncoder = commandBuffer->blitCommandEncoder();
    
    blitCommandEncoder->copyFromBuffer(sharedBuffer, 0, privateBuffer, 0, bufferLength);
    blitCommandEncoder->endEncoding();
    
    commandBuffer->commit();
    commandBuffer->waitUntilCompleted();
}

void copyDataToTexture(MTL::CommandQueue* commandQueue, MTL::Texture *privateTexture, MTL::Texture *sharedTexture) {
    MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer();
    
    MTL::BlitCommandEncoder* blitCommandEncoder = commandBuffer->blitCommandEncoder();
    
    blitCommandEncoder->copyFromTexture(sharedTexture, privateTexture);
    blitCommandEncoder->endEncoding();
    
    commandBuffer->commit();
    commandBuffer->waitUntilCompleted();
}

