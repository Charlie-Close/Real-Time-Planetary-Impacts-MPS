//
//  MatrixMath.cpp
//  SPH
//
//  Created by Charlie Close on 03/03/2025.
//

#include <stdio.h>
#include "MatrixMath.h"

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

    simd::float4x4 makeOthProj(float left, float right, float bottom, float top, float zNear, float zFar) {
        simd::float4x4 Result(1);
        Result.columns[0][0] = 2 / (right - left);
        Result.columns[1][1] = 2 / (top - bottom);
        Result.columns[2][2] = - 2 / (zFar - zNear);
        Result.columns[3][0] = - (right + left) / (right - left);
        Result.columns[3][1] = - (top + bottom) / (top - bottom);
        Result.columns[3][2] = - (zFar + zNear) / (zFar - zNear);
        return Result;
    }
}

