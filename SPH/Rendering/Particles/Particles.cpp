//
//  Particles.cpp
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

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

Particles::Particles(MTL::Device* device, MTL::Buffer* particlePositions, MTL::Buffer* extraBuffer, MTL::Buffer* densityBuffer, int nParticles) {
    nPoints = nParticles;

    buildSphereVertexBuffer(device);
    buildShaders(device);
    buildDepthStencilStates(device);

    _positionBuffer = particlePositions;
    _materialIdBuffer = extraBuffer;
    _densityBuffer = densityBuffer;
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

void Particles::draw(MTL::RenderCommandEncoder *pEnc, MTL::Buffer* cameraDataBuffer) {
    pEnc->setRenderPipelineState(_drawPSO);
    pEnc->setDepthStencilState(_depthStencilState);

    // Set buffers for instancing
    pEnc->setVertexBuffer(_sphereVertexBuffer, 0, 0); // Sphere vertex buffer
    pEnc->setVertexBuffer(_normalBuffer, 0, 1); // Normal Buffer
    pEnc->setVertexBuffer(_positionBuffer, 0, 2);    // Per-instance particle positions
    pEnc->setVertexBuffer(cameraDataBuffer, 0, 3);   // Camera data buffer
    pEnc->setVertexBuffer(_materialIdBuffer, 0, 4);       // Extra data buffer
    pEnc->setVertexBuffer(_densityBuffer, 0, 5);       // density data buffer

    // Draw instanced spheres
    pEnc->drawIndexedPrimitives( MTL::PrimitiveType::PrimitiveTypeTriangle,
                                nSphereIndices, MTL::IndexType::IndexTypeUInt16,
                                _indexBuffer,
                                0,
                                nPoints );
}
