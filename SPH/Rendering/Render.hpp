//
//  Render.hpp
//  Lattice Boltzman
//
//  Created by Charlie Close on 08/12/2024.
//

#ifndef RENDER_HPP
#define RENDER_HPP

#include <Metal/Metal.hpp>
#include <AppKit/AppKit.hpp>
#include <MetalKit/MetalKit.hpp>
#include "Camera.hpp"
#include "Particles.hpp"
#include "Compute.hpp"

class Renderer
{
public:
    Renderer( MTL::Device* pDevice, Camera* camera );
    ~Renderer();
    
    void draw( MTK::View* pView );

private:
    Camera* camera;
    Particles* particles;
    Compute* compute;
    
    
    void buildBuffers();
    
    MTL::Device* _pDevice;
    MTL::CommandQueue* _pCommandQueue;
    
    MTL::Buffer* _pCameraDataBuffer;
    
    int _frame;
    dispatch_semaphore_t _semaphore;
};


#endif
