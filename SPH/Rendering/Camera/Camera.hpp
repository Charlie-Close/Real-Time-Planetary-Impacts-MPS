//
//  Camera.hpp
//  Lattice Boltzman
//
//  Created by Charlie Close on 08/12/2024.
//

#ifndef Camera_hpp
#define Camera_hpp

#include <simd/simd.h>
#include "../../Parameters.h"

class Camera
{
    public:
        Camera();
        void handleKeyDown(int key);
        void handleKeyUp(int key);
        void handleMouseDown(float x, float y);
        void handleMouseDrag(float x, float y);
    
        simd::float4x4 getMatrix();
    

    private:
        simd_float2 lastMousePos;
        void move();
        bool downKeys[6] = { false, false, false, false, false, false };
    
        simd::float3 _position = { STARTING_POSITION };
        simd::float3 _forward = { 0.f, 0.f, 1.f };
        simd::float3 _up = { 0.f, 1.f, 0.f};
        float yaw = STARTING_YAW;
        float pitch = STARTING_PITCH;
        float near_plane = 10.f;
        float far_plane = 1000.0f;
    
        simd::float4x4 _perspective;
        simd::float4x4 _view;
};

#endif /* Camera_hpp */
