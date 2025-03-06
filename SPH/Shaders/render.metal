//
//  render.metal
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#include <metal_stdlib>
#include "../Parameters.h"
#include "utils/cellsToScan.h"
#include "utils/kernels.h"
using namespace metal;

struct VertOut {
    float4 position [[position]];
    half3 colour;
};

float calculateShadow(float4 positionInLightSpace, texture2d<float, access::sample> shadowMap) {
    positionInLightSpace.xyz /= positionInLightSpace.w;
    float2 lightSpaceCoord = positionInLightSpace.xy * 0.5 + 0.5;
    lightSpaceCoord.y = 1.0 - lightSpaceCoord.y;
    if (lightSpaceCoord.x < 0.0 || lightSpaceCoord.y < 0.0 || lightSpaceCoord.x > 1.0 || lightSpaceCoord.y > 1.0) {
        return 0;
    }
    
    sampler textureSampler(mag_filter::linear, min_filter::linear, mip_filter::none, address::clamp_to_edge);
    float lightDepth = shadowMap.sample(textureSampler, lightSpaceCoord).x;
    if (positionInLightSpace.z > lightDepth + 0.005) {
        return 1.0;
    }
    return 0;
}

vertex VertOut vertexSphere(device const float3* vertexData [[buffer(0)]],
                            device const float3* normalData [[buffer(1)]],
                            device const float3* positions [[buffer(2)]],
                            device const float3* densityGradients [[buffer(3)]],
                            constant float4x4& cameraMatrix [[buffer(4)]],
                            constant float3& cameraPos [[buffer(5)]],
                            constant float4x4& lightMatrix [[buffer(6)]],
                            device const int* materialIds [[buffer(7)]],
                            device const float* densities [[buffer(8)]],
                            device const uint* instanceIds [[buffer(9)]],
                            texture2d<float, access::sample> shadowMap [[texture(0)]],
                            uint vertexId [[vertex_id]],
                            uint instInd [[instance_id]])
{
    uint instanceId = instanceIds[instInd];
    VertOut o;
    float3 lightDir = { LIGHT_DIRECTION };
    lightDir = normalize(lightDir);

    float4 worldPosition = float4(vertexData[vertexId] + positions[instanceId], 1.0);
    o.position = cameraMatrix * worldPosition;
    float4 lightSpacePos = lightMatrix * worldPosition;
    float shadow = calculateShadow(lightSpacePos, shadowMap);
    
    if (materialIds[instanceId] == 402 or materialIds[instanceId] == 100) {
        o.colour = half3(0.2f, 0.2f, 0.25f);
    } else {
        o.colour = half3(0.5f, 0.5f, 0.6f);
    }
    
    float3 gradRho = densityGradients[instanceId];
    float gradMag = length(gradRho);
    float weight = gradMag > 1e-12 ? clamp(1000 * gradMag, 0.f, 1.f) : 0;
    float3 pNorm = - normalize(densityGradients[instanceId]);
    float3 norm = normalData[vertexId];
    float3 normal = weight == 0 ? norm : mix(norm, pNorm, weight);

    float3 viewDir = normalize(worldPosition.xyz - cameraPos);
    float3 reflectDir = reflect(-lightDir, normal);
    float spec = 0.6f * pow(max(dot(viewDir, reflectDir), 0.0), 16);
    
    o.colour *= clamp(dot(lightDir, normal), 0.f, 1.f);
    o.colour.r = min(o.colour.r + spec, 1.f);
    o.colour.g = min(o.colour.g + spec, 1.f);
    o.colour.b = min(o.colour.b + spec, 1.f);
    o.colour *= (1 - shadow);
    
    o.colour.r = clamp(o.colour.r, half(0.01), half(1));
    o.colour.g = clamp(o.colour.g, half(0.01), half(1));
    o.colour.b = clamp(o.colour.b, half(0.01), half(1));


    return o;
}

half4 fragment fragmentSphere(VertOut in [[stage_in]])
{
    return half4( in.colour, 1.0 );
}

vertex VertOut vertexShadow(device const float3* vertexData [[buffer(0)]],
                            device const float3* positions [[buffer(1)]],
                            constant float4x4& lightMatrix [[buffer(2)]],
                            device const uint* instanceIds [[buffer(3)]],
                            uint vertexId [[vertex_id]],
                            uint instInd [[instance_id]])
{
    uint instanceId = instanceIds[instInd];
    VertOut o;
    float4 worldPosition = float4(vertexData[vertexId] + positions[instanceId], 1.0);
    o.position = lightMatrix * worldPosition;
    return o;
}

half4 fragment fragmentShadow(VertOut in [[stage_in]])
{
    return half4( in.colour, 1.0 );
}

kernel void zeroVisibleCount(device uint& visibleCount,
                             uint index [[thread_position_in_grid]]) {
    visibleCount = 0;
}

kernel void setIndrBuffer(device uint& visibleCount,
                          device MTLDrawPrimitivesIndirectArguments& args,
                          uint index [[thread_position_in_grid]]) {
    args.instanceCount = visibleCount;
}

kernel void determineVisibility(device float3* positions,
                                device float* masses,
                                device float* h,
                                device float3* normals,
                                device int* cellStarts,
                                device int* cellEnds,
                                device uint2* cellData,
                                device float* densities,
                                device atomic_uint* visibleCount,
                                device uint* instanceArray,
                                constant float3& cameraPos,
                                uint index [[thread_position_in_grid]]) {
    float3 x_i = positions[index];
    float3 densityGradient = normals[index];
    float rho = densities[index];
    
    float g2 = length_squared(densityGradient);
    bool visible;
    if (g2 > 1e-12 and rho > 0.0005) {
        float g = sqrt(g2);
        float g1 = 1 / g;
        float3 normal = - densityGradient * g1;
        float distToEdge = 0.3 * (rho - 0.001) * g1;
        visible = distToEdge < 2 && (dot(x_i - cameraPos, normal) > 0 or dot(float3(LIGHT_DIRECTION), normal) > 0);
    } else {
        visible = true;
    }
    
    if (visible) {
        uint ind = atomic_fetch_add_explicit(visibleCount, 1, memory_order_relaxed);
        instanceArray[ind] = index;
    }
}
