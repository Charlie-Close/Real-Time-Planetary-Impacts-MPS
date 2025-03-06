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
    particles = new Particles(pDevice, compute->positionBuffer, compute->materialIdBuffer, compute->densityBuffer, compute->_massBuffer, compute->_smoothingLengthBuffer, compute->rhoGrads, compute->_cellStart, compute->_cellEnd, compute->_cellArrayi, compute->nParticles);

    buildBuffers();
    
    this->camera = camera;
    _semaphore = dispatch_semaphore_create(1);
}

Renderer::~Renderer()
{
    _pCameraDataBuffer->release();
    _pCommandQueue->release();
    _pDevice->release();
    delete compute;
}

void Renderer::buildBuffers()
{
    _pCameraDataBuffer = _pDevice->newBuffer( sizeof(simd::float4x4), MTL::ResourceStorageModeShared );
    _cameraPosBuffer = _pDevice->newBuffer( sizeof(simd::float3), MTL::ResourceStorageModeShared );
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
    

    // The simulation logic is in here. Can run multiple simulation steps per frame.
    for (int i = 0; i < STEPS_PER_FRAME; i++) {
        compute->sort(pCmd);
        compute->gravitationalPass(pCmd);
        if (_frame == 0) {
            // First frame loop density to get smooth gradients
            for (int j = 0; j < DENSITY_GRADIENT_SETTLING_ITTERATIONS; j++) {
                compute->densityPass(pCmd);
            }
        }
        compute->densityPass(pCmd);
        compute->accelerationPass(pCmd);
        compute->stepPass(pCmd);
        _frame++;
    }
    
    
    // Update camera matrix
    simd::float4x4* pCameraData = reinterpret_cast<simd::float4x4*>(_pCameraDataBuffer->contents());
    *pCameraData = camera->getMatrix();
    simd::float3* cameraPosition = reinterpret_cast<simd::float3*>(_cameraPosBuffer->contents());
    *cameraPosition = camera->getPosition();
    
    particles->zeroVis(pCmd);
    particles->calculateNormals(pCmd, _cameraPosBuffer);
    particles->setIndrBuff(pCmd);
    particles->drawShadowMap(pCmd);
    
    // Begin render pass
    MTL::RenderPassDescriptor* pRpd = pView->currentRenderPassDescriptor();
    MTL::RenderCommandEncoder* pEnc = pCmd->renderCommandEncoder(pRpd);
    particles->draw(pEnc, _pCameraDataBuffer, _cameraPosBuffer);
    pEnc->endEncoding();
    
    pCmd->presentDrawable(pView->currentDrawable());
    pCmd->commit();
    compute->organisation(_pDevice, _pCommandQueue);
    particles->updateBuffers(compute->positionBuffer, compute->materialIdBuffer, compute->densityBuffer, compute->_massBuffer, compute->_smoothingLengthBuffer, compute->rhoGrads);
    pPool->release();

}

