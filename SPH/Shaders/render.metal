//
//  render.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include <metal_stdlib>
#include "../Parameters.h"
using namespace metal;

struct VertOut {
    float4 position [[position]];
    half3 colour;
};

vertex VertOut vertexSphere(device const float3* vertexData [[buffer(0)]],
                            device const float3* normalData [[buffer(1)]],
                            device const float3* positions [[buffer(2)]],
                            constant float4x4& cameraMatrix [[buffer(3)]],
                            device const int* materialIds [[buffer(4)]],
                            device const float* densities [[buffer(5)]],
                            uint vertexId [[vertex_id]],
                            uint instanceId [[instance_id]])
{
    VertOut o;
    if (densities[instanceId] > DENSITY_CUTOFF) {
        o.position = { 0, 0, 0 };
        return o;
    }
    float3 lightDir = { 1.0, 0, 0 };
    float4 worldPosition = float4(vertexData[vertexId] + positions[instanceId], 1.0);
    o.position = cameraMatrix * worldPosition;
    
    if (materialIds[instanceId] == 402 or materialIds[instanceId] == 100) {
        o.colour = half3(0.2f, 0.2f, 0.25f) * max(dot(lightDir, normalData[vertexId]), 0.1);
    } else {
        o.colour = half3(0.5f, 0.5f, 0.6f) * max(dot(lightDir, normalData[vertexId]), 0.1);
    }

    return o;
}

half4 fragment fragmentSphere(VertOut in [[stage_in]])
{
    return half4( in.colour, 1.0 );
}






struct PointVertOut {
    float4 position [[position]];
    float pointSize [[point_size]];
    half3 colour;
};


vertex PointVertOut vertexPoint(device const float3* positions [[buffer(0)]],
                       constant float4x4& cameraMatrix [[buffer(1)]],
                        device const int* materialIds [[buffer(2)]],
                       uint vertexId [[vertex_id]]) {
    PointVertOut o;
    float4 worldPosition = float4(positions[vertexId], 1.0);
    o.position = cameraMatrix * worldPosition;
    
    if (materialIds[vertexId] == 400 or materialIds[vertexId] == 100) {
        o.colour = half3(1.f, 0.f, 0.f);
    } else {
        o.colour = half3(0.f, 1.f, 0.f);
    }

    // Set a defined point size (in pixels)
    o.pointSize = 1.0; // or any desired size

    return o;
}

half4 fragment fragmentPoint(PointVertOut in [[stage_in]])
{
    return half4( in.colour, 1.0 );
}

