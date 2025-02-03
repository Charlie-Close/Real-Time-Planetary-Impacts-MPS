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
    simd::float4x4 makePerspective( float fovRadians, float aspect, float znear, float zfar )
    {
        using simd::float4;
        float ys = 1.f / tanf(fovRadians * 0.5f);
        float xs = ys / aspect;
        float zs = zfar / ( znear - zfar );
        return simd_matrix_from_rows((float4){ xs, 0.0f, 0.0f, 0.0f },
                                     (float4){ 0.0f, ys, 0.0f, 0.0f },
                                     (float4){ 0.0f, 0.0f, zs, znear * zs },
                                     (float4){ 0, 0, -1, 0 });
    }
    simd::float4x4 makeLookAt(simd::float3 pos, simd::float3 forward, simd::float3 up) {
        // Normalize the forward vector
        simd::float3 zAxis = - simd::normalize(forward);

        // Compute the right vector
        simd::float3 xAxis = simd::normalize(simd::cross(up, zAxis));

        // Recompute the up vector to ensure orthogonality
        simd::float3 yAxis = simd::cross(zAxis, xAxis);

        // Create the rotation matrix
        simd::float4x4 rotation = {
            (simd::float4){xAxis.x, yAxis.x, zAxis.x, 0.0f},
            (simd::float4){xAxis.y, yAxis.y, zAxis.y, 0.0f},
            (simd::float4){xAxis.z, yAxis.z, zAxis.z, 0.0f},
            (simd::float4){0.0f, 0.0f, 0.0f, 1.0f}
        };

        // Create the translation matrix
        simd::float4x4 translation = {
            (simd::float4){1.0f, 0.0f, 0.0f, 0.0f},
            (simd::float4){0.0f, 1.0f, 0.0f, 0.0f},
            (simd::float4){0.0f, 0.0f, 1.0f, 0.0f},
            (simd::float4){-pos.x, -pos.y, -pos.z, 1.0f}
        };

        // Combine the rotation and translation matrices
        return rotation * translation;
    }
}

#endif /* MatrixMath_h */
