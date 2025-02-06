//
//  Render.cpp
//  Lattice Boltzman
//
//  Created by Charlie Close on 08/12/2024.
//

#include "Render.hpp"
#include "ViewAdapter.hpp"
#include "Buffers.hpp"
#include <iostream>

// -------------------------------- //
//                                  //
//          Constructor             //
//                                  //
// -------------------------------- //

Renderer::Renderer( MTL::Device* pDevice, Camera* camera )
: _pDevice( pDevice->retain() )
, _frame( 0 )
{
    _pCommandQueue = _pDevice->newCommandQueue();
    
    compute = new Compute(_pDevice);
    particles = new Particles(_pDevice, compute->positionBuffer, compute->_materialIdBuffer, compute->nParticles);

    buildBuffers();
    
    this->camera = camera;
    _semaphore = dispatch_semaphore_create(1);
}

Renderer::~Renderer()
{
    _pCameraDataBuffer->release();
    _pCommandQueue->release();
    _pDevice->release();
}

void Renderer::buildBuffers()
{
    _pCameraDataBuffer = _pDevice->newBuffer( sizeof(simd::float4x4), MTL::ResourceStorageModeManaged );
}

// -------------------------------- //
//                                  //
//          Draw command            //
//                                  //
// -------------------------------- //

void Renderer::draw(MTK::View* pView)
{
    NS::AutoreleasePool* pPool = NS::AutoreleasePool::alloc()->init();
    MTL::CommandBuffer* pCmd = _pCommandQueue->commandBuffer();
    
    // Semaphore stuff
    dispatch_semaphore_wait(_semaphore, DISPATCH_TIME_FOREVER);
    Renderer* pRenderer = this;
    pCmd->addCompletedHandler(^void(MTL::CommandBuffer* pCmd) {
        dispatch_semaphore_signal(pRenderer->_semaphore);
    });
    

    compute->updateOctreeBuffer(_pDevice);
    for (int i = 0; i < 2; i++) {
        compute->sort(pCmd);
        compute->densityPass(pCmd);
        compute->gravitationalPass(pCmd);
        compute->accelerationPass(pCmd);
        compute->stepPass(pCmd);
        _frame++;
    }
    
    // Update camera matrix
    simd::float4x4* pCameraData = reinterpret_cast<simd::float4x4*>(_pCameraDataBuffer->contents());
    *pCameraData = camera->getMatrix();
    _pCameraDataBuffer->didModifyRange(NS::Range::Make(0, sizeof(simd::float4x4)));
    
    // Begin render pass
    MTL::RenderPassDescriptor* pRpd = pView->currentRenderPassDescriptor();
    MTL::RenderCommandEncoder* pEnc = pCmd->renderCommandEncoder(pRpd);
    
    particles->draw(pEnc, _pCameraDataBuffer);

    pEnc->endEncoding();
    
    pCmd->presentDrawable(pView->currentDrawable());
    pCmd->commit();

    pPool->release();
}

