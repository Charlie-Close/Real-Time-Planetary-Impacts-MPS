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
#include "Snapshotter.hpp"
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
    
    int nextSnapshot = 0;
    int _frame = 0;
    
    Compute* compute;
    Particles* particles;
    Snapshotter** snapshotters;
    Camera* camera;
    
    MTL::Device* _device;
    MTL::CommandQueue* _commandQueue;
    
    MTL::Buffer* _cameraDataBuffer;
    MTL::Buffer* _cameraPosBuffer;
};

#endif /* Headless_hpp */
