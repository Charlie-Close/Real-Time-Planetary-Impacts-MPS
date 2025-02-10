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
    // Cull vertices that are highly dense (assume they are not on surface, and therefore not visible)
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

kernel void drawToSnapshot(device float3* positions,
                           device int* materialIds,
                           device atomic_int* output,
                           uint index [[thread_position_in_grid]])
{
    const float R_E = 6.3710;
    float2 uv = positions[index].xy;// / (20 * R_E);
    uv.x -= 250;
    uv.y -= 250;
    uv *= SNAPSHOT_RESOLUTION / (20 * R_E);
    if (uv.x < 0 or uv.y < 0 or uv.x > SNAPSHOT_RESOLUTION or uv.y > SNAPSHOT_RESOLUTION) {
        return;
    }
    
    int pixel = (int)uv.x + (SNAPSHOT_RESOLUTION - 1 - (int)uv.y) * SNAPSHOT_RESOLUTION;
    
    atomic_fetch_max_explicit(&(output[pixel]), materialIds[index], memory_order_relaxed);
}
