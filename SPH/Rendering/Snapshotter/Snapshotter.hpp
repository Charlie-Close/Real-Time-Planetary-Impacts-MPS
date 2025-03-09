//
//  Snapshotter.hpp
//  SPH
//
//  Created by Charlie Close on 09/03/2025.
//

#ifndef Snapshotter_hpp
#define Snapshotter_hpp

#include <stdio.h>
#include <Metal/Metal.hpp>
#include "Compute.hpp"
#include "Particles.hpp"
#include "Camera.hpp"
#include <vector>
#include <string>

class Snapshotter
{
public:
    Snapshotter(MTL::Device* device, Particles* particles, simd_float3 position, float pitch, float yaw, std::string snapshotFolder);
    ~Snapshotter();
    void takeSnapshot(int n);

private:
    bool loading = false;
    void createRenderPassDescriptor();
    std::string _snapshotFolder;
    
    int nextSnapshot = 0;
    int _frame = 0;
    
    Particles* _particles;
    Camera* _camera;
    
    MTL::Device* _device;
    MTL::CommandQueue* _commandQueue;
    
    MTL::Texture* colourTexture;
    MTL::Texture* depthTexture;
    MTL::Buffer* _cameraDataBuffer;
    MTL::Buffer* _cameraPosBuffer;
    
    MTL::RenderPassDescriptor* renderPassDescriptor;
};


#endif /* Snapshotter_hpp */
