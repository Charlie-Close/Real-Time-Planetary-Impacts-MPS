//
//  Headless.hpp
//  SPH
//
//  Created by Charlie Close on 10/02/2025.
//

#ifndef Headless_hpp
#define Headless_hpp

#include <stdio.h>
#include <Metal/Metal.hpp>
#include "Compute.hpp"
#include "Particles.hpp"
#include "Camera.hpp"
#include <vector>

class Headless
{
public:
    Headless();
    ~Headless();
    
    void step();
    float getTime();

private:
    bool takeSnapshot();
    bool loading = false;
    void createRenderPassDescriptor();
    std::pair<MTL::Texture*, MTL::Texture*> createTextures();
    
    int nextSnapshot = 0;
    int _frame = 0;
    
    Compute* compute;
    Particles* particles;
    Camera* camera;
    
    MTL::Device* _device;
    MTL::CommandQueue* _commandQueue;
    
    MTL::Texture* colourTexture;
    MTL::Texture* depthTexture;
    MTL::Buffer* _cameraDataBuffer;
    MTL::Buffer* _cameraPosBuffer;
    
    MTL::RenderPassDescriptor* renderPassDescriptor;
};

#endif /* Headless_hpp */
