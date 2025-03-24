//
//  Particles.cpp
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include "MatrixMath.h"
#include "Particles.hpp"
#include "Buffers.hpp"
#include "ParticleMesh.hpp"
#include "../../Parameters.h"
#include <iostream>

// -------------------------------- //
//                                  //
//          Constructor             //
//                                  //
// -------------------------------- //

Particles::Particles(MTL::Device* device, MTL::Buffer* particlePositions, MTL::Buffer* materialIdBuffer, MTL::Buffer* densityBuffer, MTL::Buffer* massBuffer, MTL::Buffer* smoothingLengthBuffer, MTL::Buffer* rhoGrads, MTL::Buffer* temperature, MTL::Buffer* alpha, int nParticles) {
    nPoints = nParticles;
    
    buildSphereVertexBuffer(device);
    buildShaders(device);
    buildDepthStencilStates(device);
    
    _positionBuffer = particlePositions;
    _materialIdBuffer = materialIdBuffer;
    _densityBuffer = densityBuffer;
    _massBuffer = massBuffer;
    _smoothingLengthBuffer = smoothingLengthBuffer;
    _pNormalBuffer = rhoGrads;
    _temperature = temperature;
    _alpha = alpha;
    
    NS::Error** error = nil;
    MTL::Library* defaultLibrary = device->newDefaultLibrary();
    MTL::Function* visFunction = defaultLibrary->newFunction(NS::String::string("determineVisibility", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* zeroFunction = defaultLibrary->newFunction(NS::String::string("zeroVisibleCount", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* setIndrFunction = defaultLibrary->newFunction(NS::String::string("setIndrBuffer", NS::StringEncoding::UTF8StringEncoding));

    _visPSO = device->newComputePipelineState(visFunction, error);
    _zeroVisiblePSO = device->newComputePipelineState(zeroFunction, error);
    _setIndrBuffPso = device->newComputePipelineState(setIndrFunction, error);

    
    MTL::TextureDescriptor* depthDesc = MTL::TextureDescriptor::texture2DDescriptor(
        MTL::PixelFormatDepth32Float,
        SHADOW_MAP_RESOLUTION,
        SHADOW_MAP_RESOLUTION,
        false
    );
    depthDesc->setUsage(MTL::TextureUsageRenderTarget);
    depthDesc->setUsage(MTL::TextureUsageShaderRead);
    _shadowMap = device->newTexture(depthDesc);
    
    shadowMapPassDescriptor = MTL::RenderPassDescriptor::renderPassDescriptor();
    MTL::RenderPassDepthAttachmentDescriptor* depthAttachementDescriptor = shadowMapPassDescriptor->depthAttachment();
    depthAttachementDescriptor->setTexture(_shadowMap);
    depthAttachementDescriptor->setLoadAction(MTL::LoadActionClear);
}

void Particles::updateBuffers(MTL::Buffer* particlePositions, MTL::Buffer* materialIdBuffer, MTL::Buffer* densityBuffer, MTL::Buffer* massBuffer, MTL::Buffer* smoothingLengthBuffer, MTL::Buffer* rhoGrads, MTL::Buffer* temperature, MTL::Buffer* alpha) {
    _positionBuffer = particlePositions;
    _materialIdBuffer = materialIdBuffer;
    _densityBuffer = densityBuffer;
    _massBuffer = massBuffer;
    _smoothingLengthBuffer = smoothingLengthBuffer;
    _pNormalBuffer = rhoGrads;
    _temperature = temperature;
    _alpha = alpha;
}

// -------------------------------- //
//                                  //
//      Metal Object Builders       //
//                                  //
// -------------------------------- //

void Particles::buildSphereVertexBuffer(MTL::Device* device) {
    std::vector<simd_float3> vertices;
    std::vector<simd_float3> normals;
    std::vector<uint16_t> indices;
    
    std::tie(vertices, normals, indices) = generateSphere(PARTICLE_SIZE, PARTICLE_SUBDIVITIONS);
    nSphereIndices = (int)indices.size() * 3;
    
    _sphereVertexBuffer = device->newBuffer( vertices.size() * sizeof(simd_float3), MTL::ResourceStorageModeShared );
    _normalBuffer = device->newBuffer( normals.size() * sizeof(simd_float3), MTL::ResourceStorageModeShared );
    _indexBuffer = device->newBuffer( indices.size() * sizeof(uint16_t) * 3, MTL::ResourceStorageModeShared );
    _indrBuffer = device->newBuffer(sizeof(MTL::DrawIndexedPrimitivesIndirectArguments), MTL::ResourceStorageModePrivate);
    _lightMatrixBuffer = device->newBuffer(sizeof(simd::float4x4), MTL::ResourceStorageModePrivate);
    _visibleCount = device->newBuffer(sizeof(uint), MTL::ResourceStorageModePrivate);
    _instanceArray = device->newBuffer(nPoints * sizeof(uint), MTL::ResourceStorageModePrivate);

    MTL::CommandQueue* cmdQueue = device->newCommandQueue();
    simd::float3 lightDirection = { LIGHT_DIRECTION };
    lightDirection = simd::normalize(lightDirection);
    simd::float3 center(BOX_CENTER);
    simd::float4x4 view = math::makeLookAt(center - lightDirection * 300, lightDirection, { 0, 1, 0 });
    simd::float4x4 proj = math::makeOthProj(-50, 50, -50, 50, 0, 500);
    simd::float4x4 lightMatrix = proj * view;
    writeDataToPrivateBuffer(device, cmdQueue, _lightMatrixBuffer, lightMatrix);
    MTL::DrawIndexedPrimitivesIndirectArguments args;
    args.indexCount = nSphereIndices;
    args.instanceCount = nPoints;
    args.indexStart = 0;
    args.baseVertex = 0;
    args.baseInstance = 0;
    writeDataToPrivateBuffer(device, cmdQueue, _indrBuffer, args);
    
    writeDataToBuffer(_sphereVertexBuffer, vertices);
    writeDataToBuffer(_normalBuffer, normals);
    writeDataToBuffer(_indexBuffer, indices);
}

void Particles::buildShaders(MTL::Device* device) {
    NS::Error** error = nil;

    MTL::Library* defaultLibrary = device->newDefaultLibrary();

    MTL::Function* pVertexFn = defaultLibrary->newFunction(NS::String::string("vertexSphere", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* pFragFn = defaultLibrary->newFunction(NS::String::string("fragmentSphere", NS::StringEncoding::UTF8StringEncoding));

    MTL::RenderPipelineDescriptor* pDesc = MTL::RenderPipelineDescriptor::alloc()->init();
    pDesc->setVertexFunction(pVertexFn);
    pDesc->setFragmentFunction(pFragFn);
    pDesc->colorAttachments()->object(0)->setPixelFormat(MTL::PixelFormat::PixelFormatBGRA8Unorm_sRGB);
    pDesc->setDepthAttachmentPixelFormat(MTL::PixelFormat::PixelFormatDepth16Unorm);

    _drawPSO = device->newRenderPipelineState(pDesc, error);

    pVertexFn->release();
    pFragFn->release();
    pDesc->release();
    
    MTL::Function* vShadFn = defaultLibrary->newFunction(NS::String::string("vertexShadow", NS::StringEncoding::UTF8StringEncoding));
    MTL::Function* fShadFn = defaultLibrary->newFunction(NS::String::string("fragmentShadow", NS::StringEncoding::UTF8StringEncoding));

    MTL::RenderPipelineDescriptor* shadDesc = MTL::RenderPipelineDescriptor::alloc()->init();
    shadDesc->setVertexFunction(vShadFn);
    shadDesc->setFragmentFunction(fShadFn);
    shadDesc->setDepthAttachmentPixelFormat(MTL::PixelFormat::PixelFormatDepth32Float);
    
    _shadowPSO = device->newRenderPipelineState(shadDesc, error);

    vShadFn->release();
    fShadFn->release();
    shadDesc->release();
}

void Particles::buildDepthStencilStates(MTL::Device* device) {
    MTL::DepthStencilDescriptor* pDsDesc = MTL::DepthStencilDescriptor::alloc()->init();
    pDsDesc->setDepthCompareFunction(MTL::CompareFunction::CompareFunctionLess);
    pDsDesc->setDepthWriteEnabled(true);
    _depthStencilState = device->newDepthStencilState(pDsDesc);
    pDsDesc->release();
}

// -------------------------------- //
//                                  //
//              Drawing             //
//                                  //
// -------------------------------- //

void Particles::zeroVis(MTL::CommandBuffer* cmdBuff) {
    MTL::ComputeCommandEncoder* cmdEncoder = cmdBuff->computeCommandEncoder();
    cmdEncoder->setBuffer(_visibleCount, 0, 0);
    MTL::Size gridSize = MTL::Size(1, 1, 1);
    cmdEncoder->setComputePipelineState(_zeroVisiblePSO);
    MTL::Size threadgroupSize = MTL::Size(1, 1, 1);
    cmdEncoder->dispatchThreads(gridSize, threadgroupSize);
    cmdEncoder->endEncoding();
}
void Particles::setIndrBuff(MTL::CommandBuffer* cmdBuff) {
    MTL::ComputeCommandEncoder* cmdEncoder = cmdBuff->computeCommandEncoder();
    cmdEncoder->setBuffer(_visibleCount, 0, 0);
    cmdEncoder->setBuffer(_indrBuffer, 0, 1);
    MTL::Size gridSize = MTL::Size(1, 1, 1);
    cmdEncoder->setComputePipelineState(_setIndrBuffPso);
    MTL::Size threadgroupSize = MTL::Size(1, 1, 1);
    cmdEncoder->dispatchThreads(gridSize, threadgroupSize);
    cmdEncoder->endEncoding();
}

void Particles::calculateNormals(MTL::CommandBuffer* cmdBuff, MTL::Buffer* cameraPosBuffer) {
    MTL::ComputeCommandEncoder* cmdEncoder = cmdBuff->computeCommandEncoder();
    cmdEncoder->setBuffer(_positionBuffer, 0, 0);
    cmdEncoder->setBuffer(_massBuffer, 0, 1);
    cmdEncoder->setBuffer(_smoothingLengthBuffer, 0, 2);
    cmdEncoder->setBuffer(_pNormalBuffer, 0, 3);
    cmdEncoder->setBuffer(_densityBuffer, 0, 4);
    cmdEncoder->setBuffer(_visibleCount, 0, 5);
    cmdEncoder->setBuffer(_instanceArray, 0, 6);
    cmdEncoder->setBuffer(cameraPosBuffer, 0, 7);
    MTL::Size gridSize = MTL::Size(nPoints, 1, 1);
    cmdEncoder->setComputePipelineState(_visPSO);
    NS::UInteger threadGroupSize = fmin(_visPSO->maxTotalThreadsPerThreadgroup(), nPoints);
    MTL::Size threadgroupSize = MTL::Size(threadGroupSize, 1, 1);
    cmdEncoder->dispatchThreads(gridSize, threadgroupSize);
    cmdEncoder->endEncoding();
}

void Particles::drawShadowMap(MTL::CommandBuffer *cmdBuffer) {
    shadowMapPassDescriptor = MTL::RenderPassDescriptor::renderPassDescriptor();
    MTL::RenderPassDepthAttachmentDescriptor* depthAttachementDescriptor = shadowMapPassDescriptor->depthAttachment();
    depthAttachementDescriptor->setTexture(_shadowMap);
    depthAttachementDescriptor->setLoadAction(MTL::LoadActionClear);
    depthAttachementDescriptor->setClearDepth(1.0);
    depthAttachementDescriptor->setStoreAction(MTL::StoreActionStore);
    MTL::RenderCommandEncoder* cmdEnc = cmdBuffer->renderCommandEncoder(shadowMapPassDescriptor);
    cmdEnc->setRenderPipelineState(_shadowPSO);
    cmdEnc->setDepthStencilState(_depthStencilState);
    cmdEnc->setVertexBuffer(_sphereVertexBuffer, 0, 0); // Sphere vertex buffer
    cmdEnc->setVertexBuffer(_positionBuffer, 0, 1);    // Per-instance particle positions
    cmdEnc->setVertexBuffer(_lightMatrixBuffer, 0, 2);   // Camera data buffer
    cmdEnc->setVertexBuffer(_instanceArray, 0, 3);
    cmdEnc->setVertexBuffer(_densityBuffer, 0, 4);
    cmdEnc->setVertexBuffer(_smoothingLengthBuffer, 0, 5);
    cmdEnc->drawIndexedPrimitives(MTL::PrimitiveType::PrimitiveTypeTriangle, MTL::IndexType::IndexTypeUInt16, _indexBuffer, 0, _indrBuffer, 0);
    cmdEnc->endEncoding();
}

void Particles::draw(MTL::RenderCommandEncoder *pEnc, MTL::Buffer* cameraDataBuffer, MTL::Buffer* cameraPosBuffer) {
    pEnc->setRenderPipelineState(_drawPSO);
    pEnc->setDepthStencilState(_depthStencilState);

    // Set buffers for instancing
    pEnc->setVertexBuffer(_sphereVertexBuffer, 0, 0); // Sphere vertex buffer
    pEnc->setVertexBuffer(_normalBuffer, 0, 1); // Normal Buffer
    pEnc->setVertexBuffer(_positionBuffer, 0, 2);    // Per-instance particle positions
    pEnc->setVertexBuffer(_pNormalBuffer, 0, 3);    // Per-instance particle normals
    pEnc->setVertexBuffer(cameraDataBuffer, 0, 4);   // Camera data buffer
    pEnc->setVertexBuffer(cameraPosBuffer, 0, 5);   // Camera data buffer
    pEnc->setVertexBuffer(_lightMatrixBuffer, 0, 6);   // Camera data buffer
    pEnc->setVertexBuffer(_materialIdBuffer, 0, 7);       // Extra data buffer
    pEnc->setVertexBuffer(_densityBuffer, 0, 8);       // density data buffer
    pEnc->setVertexBuffer(_instanceArray, 0, 9);
    pEnc->setVertexBuffer(_temperature, 0, 10);
    pEnc->setVertexBuffer(_smoothingLengthBuffer, 0, 11);
    pEnc->setVertexBuffer(_alpha, 0, 12);
    pEnc->setFragmentTexture(_shadowMap, 0);

    pEnc->drawIndexedPrimitives(MTL::PrimitiveType::PrimitiveTypeTriangle, MTL::IndexType::IndexTypeUInt16, _indexBuffer, 0, _indrBuffer, 0);
}
