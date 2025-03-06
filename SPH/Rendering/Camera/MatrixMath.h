//
//  MatrixMath.h
//  Lattice Boltzman
//
//  Created by Charlie Close on 19/01/2025.
//

#ifndef MatrixMath_h
#define MatrixMath_h

#include <simd/simd.h>

namespace math
{
    simd::float4x4 makePerspective(float fovRadians, float aspect, float znear, float zfar);
    simd::float4x4 makeLookAt(simd::float3 pos, simd::float3 forward, simd::float3 up);
    simd::float4x4 makeOthProj(float left, float right, float bottom, float top, float zNear, float zFar);
}

#endif /* MatrixMath_h */
