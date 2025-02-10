//
//  Headless.cpp
//  SPH
//
//  Created by Charlie Close on 10/02/2025.
//

#include "Headless.hpp"
#include "Buffers.hpp"
#include <iostream>
#include <chrono>

// -------------------------------- //
//                                  //
//          Constructor             //
//                                  //
// -------------------------------- //

Headless::Headless() {
    _device = MTL::CreateSystemDefaultDevice();
    _commandQueue = _device->newCommandQueue();
    compute = new Compute(_device);
}

Headless::~Headless()
{
    _commandQueue->release();
    _device->release();
}

// -------------------------------- //
//                                  //
//          Draw command            //
//                                  //
// -------------------------------- //

void Headless::step()
{
    MTL::CommandBuffer* pCmd = _commandQueue->commandBuffer();
    
    // The simulation logic is in here. Can run multiple simulation steps per frame.
    // Only attempt to update octree once (this is running on CPU, so running multiple times won't do anything.
    compute->updateOctreeBuffer(_device);
    compute->drawSnapshot(pCmd);
    for (int i = 0; i < STEPS_PER_FRAME; i++) {
        compute->sort(pCmd);
        compute->gravitationalPass(pCmd);
        compute->densityPass(pCmd);
        compute->accelerationPass(pCmd);
        compute->stepPass(pCmd);
    }
    
    pCmd->commit();
    pCmd->waitUntilCompleted();
}


