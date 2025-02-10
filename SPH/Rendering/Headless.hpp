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

class Headless
{
public:
    Headless();
    ~Headless();
    
    void step();

private:
    Compute* compute;
    
    MTL::Device* _device;
    MTL::CommandQueue* _commandQueue;
};

#endif /* Headless_hpp */
