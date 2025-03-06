//
//  Particles.hpp
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#ifndef Particles_hpp
#define Particles_hpp

#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <QuartzCore/QuartzCore.hpp>
#include <simd/simd.h>
#include <string>
#include <vector>
#include <AppKit/AppKit.hpp>
#include <MetalKit/MetalKit.hpp>

class Particles {
public:
    Particles(MTL::Device* _device, MTL::Buffer* particlePositions, MTL::Buffer* materialIdBuffer, MTL::Buffer* densityBuffer, MTL::Buffer* _massBuffer, MTL::Buffer* _smoothingLengthBuffer, MTL::Buffer* rhoGrads, MTL::Buffer* _cellStarts, MTL::Buffer* _cellEnds, MTL::Buffer* _cellData, int nParticles);
    void zeroVis(MTL::CommandBuffer* cmdBuff);
    void calculateNormals(MTL::CommandBuffer* cmdBuff, MTL::Buffer* cameraPosBuffer);
    void setIndrBuff(MTL::CommandBuffer* cmdBuff);
    void drawShadowMap(MTL::CommandBuffer* cmdBuffer);
    void draw(MTL::RenderCommandEncoder *pEnc, MTL::Buffer* cameraDataBuffer, MTL::Buffer* cameraPosBuffer);
    void updateBuffers(MTL::Buffer* particlePositions, MTL::Buffer* materialIdBuffer, MTL::Buffer* densityBuffer, MTL::Buffer* massBuffer, MTL::Buffer* smoothingLengthBuffer, MTL::Buffer* rhoGrads);


private:
    void buildSphereVertexBuffer(MTL::Device* device);
    void buildShaders(MTL::Device* device);
    void buildDepthStencilStates(MTL::Device* device);
    void encodeCommand(MTL::ComputeCommandEncoder* computeEncoder, MTL::ComputePipelineState* command);
    
    MTL::Buffer* _sphereVertexBuffer;
    MTL::Buffer* _normalBuffer;
    MTL::Buffer* _indexBuffer;
    MTL::Buffer* _pNormalBuffer;
    
    int nSphereIndices;
    int nPoints;

    // Pipleline State Objects
    MTL::RenderPipelineState* _shadowPSO;
    MTL::RenderPipelineState* _drawPSO;
    MTL::ComputePipelineState* _visPSO;
    MTL::ComputePipelineState* _zeroVisiblePSO;
    MTL::ComputePipelineState* _setIndrBuffPso;

    
    // Depth stencil
    MTL::DepthStencilState* _depthStencilState;
    
    // Buffers
    MTL::Buffer* _positionBuffer;
    MTL::Buffer* _materialIdBuffer;
    MTL::Buffer* _densityBuffer;
    MTL::Buffer* _massBuffer;
    MTL::Buffer* _smoothingLengthBuffer;
    MTL::Buffer* _cellStarts;
    MTL::Buffer* _cellEnds;
    MTL::Buffer* _cellData;
    MTL::Buffer* _lightMatrixBuffer;
    MTL::Buffer* _indrBuffer;
    MTL::Buffer* _visibleCount;
    MTL::Buffer* _instanceArray;
    
    // Textures
    MTL::Texture* _shadowMap;
    
    MTL::RenderPassDescriptor* shadowMapPassDescriptor;
};

#endif /* Particles_hpp */
