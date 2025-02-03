//
//  ViewDelegate.hpp
//  Lattice Boltzman
//
//  Created by Charlie Close on 11/12/2024.
//

#ifndef ViewDelegate_hpp
#define ViewDelegate_hpp

#include "Camera.hpp"
#include "Render.hpp"

class MyMTKViewDelegate : public MTK::ViewDelegate
{
    public:
        MyMTKViewDelegate( MTL::Device* pDevice, Camera* camera );
        virtual ~MyMTKViewDelegate() override;
        virtual void drawInMTKView( MTK::View* pView ) override;
        Renderer* _pRenderer;
};

#endif /* ViewDelegate_hpp */
