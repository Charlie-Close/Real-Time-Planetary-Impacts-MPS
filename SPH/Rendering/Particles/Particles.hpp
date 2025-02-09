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
    Particles(MTL::Device* _device, MTL::Buffer* particlePositions, MTL::Buffer* extraBuffer, MTL::Buffer* densityBuffer, int nParticles);
    void draw(MTL::RenderCommandEncoder *pEnc, MTL::Buffer* cameraDataBuffer);


private:
    void buildSphereVertexBuffer(MTL::Device* device);
    void buildShaders(MTL::Device* device);
    void buildDepthStencilStates(MTL::Device* device);
    void encodeCommand(MTL::ComputeCommandEncoder* computeEncoder, MTL::ComputePipelineState* command);
    
    MTL::Buffer* _sphereVertexBuffer;
    MTL::Buffer* _normalBuffer;
    MTL::Buffer* _indexBuffer;
    
    int nSphereIndices;
    int nPoints;

    // Pipleline State Objects
    MTL::RenderPipelineState* _drawPSO;
    
    // Depth stencil
    MTL::DepthStencilState* _depthStencilState;
    
    // Buffers
    MTL::Buffer* _positionBuffer;
    MTL::Buffer* _materialIdBuffer;
    MTL::Buffer* _densityBuffer;
};

#endif /* Particles_hpp */
