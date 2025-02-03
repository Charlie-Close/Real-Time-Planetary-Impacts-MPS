//
//  ViewDelegate.cpp
//  Lattice Boltzman
//
//  Created by Charlie Close on 11/12/2024.
//

#include "ViewDelegate.hpp"

MyMTKViewDelegate::MyMTKViewDelegate( MTL::Device* pDevice, Camera* camera )
: MTK::ViewDelegate()
, _pRenderer( new Renderer( pDevice, camera ) )
{
}

MyMTKViewDelegate::~MyMTKViewDelegate()
{
    delete _pRenderer;
}

void MyMTKViewDelegate::drawInMTKView( MTK::View* pView )
{
    _pRenderer->draw( pView );
}

