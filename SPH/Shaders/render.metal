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
    float2 vertPos;
    float3 right;
    float3 up;
    float3 viewNorm;
    float3 rhoGradNorm;
    float4 lightspacePos;
    float weight;
    half3 colour;
    half3 blackbody;
};

half3 getBlackBody(float T) {
    // T is temperature in Kelvin
    float t = T / 1000.0f; // scale temperature
    float r, g, b;

    if (T < 1000.0f) T = 1000.0f;    // clamp to avoid out-of-range
    if (T > 40000.0f) T = 40000.0f;

    // Red
    if (t <= 66.0f) {
        r = 255.0f;
    } else {
        float tmp = t - 60.0f;
        tmp = 329.698727446f * pow(tmp, -0.1332047592f);
        r = clamp(tmp, 0.0f, 255.0f);
    }

    // Green
    if (t <= 66.0f) {
        float tmp = t;
        tmp = 99.4708025861f * log(tmp) - 161.1195681661f;
        g = clamp(tmp, 0.0f, 255.0f);
    } else {
        float tmp = t - 60.0f;
        tmp = 288.1221695283f * pow(tmp, -0.0755148492f);
        g = clamp(tmp, 0.0f, 255.0f);
    }

    // Blue
    if (t >= 66.0f) {
        b = 255.0f;
    } else {
        if (t <= 19.0f) {
            b = 0.0f;
        } else {
            float tmp = t - 10.0f;
            tmp = 138.5177312231f * log(tmp) - 305.0447927307f;
            b = clamp(tmp, 0.0f, 255.0f);
        }
    }

    // Convert from 0-255 range to 0-1 range
    return half3(r, g, b) / 255.0;
}

float calculateShadow(float4 positionInLightSpace, texture2d<float, access::sample> shadowMap) {
    positionInLightSpace.xyz /= positionInLightSpace.w;
    float2 lightSpaceCoord = positionInLightSpace.xy * 0.5 + 0.5;
    lightSpaceCoord.y = 1.0 - lightSpaceCoord.y;
    if (lightSpaceCoord.x < 0.0 || lightSpaceCoord.y < 0.0 || lightSpaceCoord.x > 1.0 || lightSpaceCoord.y > 1.0) {
        return 0;
    }
    
    sampler textureSampler(mag_filter::linear, min_filter::linear, mip_filter::none, address::clamp_to_edge);
    float lightDepth = shadowMap.sample(textureSampler, lightSpaceCoord).x;
    if (positionInLightSpace.z > lightDepth) {
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
                            device const float* temperatures [[buffer(10)]],
                            device const float* h [[buffer(11)]],
                            device const float* alpha [[buffer(12)]],
                            uint vertexId [[vertex_id]],
                            uint instInd [[instance_id]])
{
    uint instanceId = instanceIds[instInd];
    VertOut o;
    o.vertPos = vertexData[vertexId].xy;
    float3 lightDir = { LIGHT_DIRECTION };
    lightDir = normalize(lightDir);
    float3 viewNorm = normalize(positions[instanceId] - cameraPos);
    float3 right = normalize(cross(viewNorm, float3(0, 1, 0)));
    float3 up = normalize(cross(viewNorm, right));
    o.right = right;
    o.up = up;
    o.viewNorm = viewNorm;
    
    
    float3 gradRho = densityGradients[instanceId];
    float gradMag = length(gradRho);
    float rho = densities[instanceId];
    float weight = gradMag > 1e-12 ? clamp(2000 * (rho - 0.001), 0.f, 1.f) : 0;
//    float weight = gradMag > 1e-12 ? clamp(2000 * rho, 0.f, 1.f) : 0;
    o.weight = min(weight, 1.f);
    o.rhoGradNorm = gradRho / gradMag;
    float size = PARTICLE_SIZE * rho * h[instanceId];
    
    float3 vertPos = size * (vertexData[vertexId].x * right + vertexData[vertexId].y * up);
    if (rho < 0.0005) vertPos.x = 100000; // Don't bother rendering under dense particles
    
    float4 worldPosition = float4(vertPos + positions[instanceId], 1.0);
    o.position = cameraMatrix * worldPosition;
    o.lightspacePos = lightMatrix * (worldPosition - 8 * float4(lightDir, 0) * size); // step in light direction to stop self shadowing
    
    switch (materialIds[instanceId]) {
        case (201): {
            o.colour = half3(0.05, 0.05, .7);
            break;
        } default: {
            o.colour = half3(0.3f, 0.3f, 0.35f);
            break;
        }
    }
    
    float T = temperatures[instanceId];
//    float T = 8000 * alpha[instanceId];
//    float T = 100000 * densities[instanceId];
    float emmision = pow(clamp(T / 8000, 0.f, 1.f), 4);
    o.blackbody = emmision * getBlackBody(T);
    
    return o;
}

half4 fragment fragmentSphere(texture2d<float, access::sample> shadowMap [[texture(0)]],
                              VertOut in [[stage_in]])
{
    float r_2 = length_squared(in.vertPos);
    if (r_2 > 1.f) {
        discard_fragment();
    } else {
        float3 lightDir = { LIGHT_DIRECTION };
        lightDir = normalize(lightDir);
        float z = sqrt(1 - r_2);
        float3 normal = in.vertPos.x * in.right + in.vertPos.y * in.up - z * in.viewNorm;
        if (in.weight != 0) normal = mix(normal, in.rhoGradNorm, in.weight);
        in.colour *= max(dot(normal, -lightDir), 0.f);
        
        float3 reflectDir = reflect(- lightDir, normal);
        float spec = 0.6f * pow(max(dot(in.viewNorm, reflectDir), 0.0), 16);
        in.colour.r = min(in.colour.r + spec, 1.f);
        in.colour.g = min(in.colour.g + spec, 1.f);
        in.colour.b = min(in.colour.b + spec, 1.f);
        
        float shadow = calculateShadow(in.lightspacePos, shadowMap);
        in.colour *= (1 - shadow);
        
        in.colour += in.blackbody;
    
        in.colour.r = clamp(in.colour.r, half(0.01), half(1));
        in.colour.g = clamp(in.colour.g, half(0.01), half(1));
        in.colour.b = clamp(in.colour.b, half(0.01), half(1));
    }
    return half4( in.colour, 1.0 );
}

vertex VertOut vertexShadow(device const float3* vertexData [[buffer(0)]],
                            device const float3* positions [[buffer(1)]],
                            constant float4x4& lightMatrix [[buffer(2)]],
                            device const uint* instanceIds [[buffer(3)]],
                            device const float* densities [[buffer(4)]],
                            device const float* h [[buffer(5)]],
                            uint vertexId [[vertex_id]],
                            uint instInd [[instance_id]])
{
    uint instanceId = instanceIds[instInd];
    VertOut o;
    o.vertPos = vertexData[vertexId].xy;
    float3 ld = float3(LIGHT_DIRECTION);
    float3 right = normalize(cross(ld, float3(0, 1, 0)));
    float3 up = normalize(cross(ld, right));
    float rho = densities[instanceId];
    float size = PARTICLE_SIZE * rho * h[instanceId];
    float3 vertPos = size * (vertexData[vertexId].x * right + vertexData[vertexId].y * up);
    if (rho < 0.0005) vertPos.x = 100000;
    
    float4 worldPosition = float4(vertPos + positions[instanceId], 1.0);
    o.position = lightMatrix * worldPosition;
    return o;
}

half4 fragment fragmentShadow(VertOut in [[stage_in]])
{
    if (length_squared(in.vertPos) > 1.f) {
        discard_fragment();
    }
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
                                device float* densities,
                                device atomic_uint* visibleCount,
                                device uint* instanceArray,
                                constant float3& cameraPos,
                                uint index [[thread_position_in_grid]]) {
    float3 x_i = positions[index];
    float h_i = h[index];
    float3 densityGradient = normals[index];
    float rho = densities[index];
    
    float g2 = length_squared(densityGradient);
//    bool visible = x_i.z > BOX_SIZE;
    bool visible = true;
    if (g2 > 1e-12 and rho > 0.0005) {
        float g = sqrt(g2);
        float g1 = 1 / g;
        float3 normal = - densityGradient * g1;
        float distToEdge = 0.3 * (rho - 0.001) * g1;
        visible = distToEdge < 4 * sqrt(h_i) && (dot(x_i - cameraPos, normal) > 0 or dot(float3(LIGHT_DIRECTION), normal) > 0);
    } else {
        visible = true;
    }
    
    if (visible) {
        uint ind = atomic_fetch_add_explicit(visibleCount, 1, memory_order_relaxed);
        instanceArray[ind] = index;
    }
}
