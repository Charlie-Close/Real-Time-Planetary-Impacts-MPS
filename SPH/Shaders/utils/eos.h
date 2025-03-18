//
//  eos.h
//  SPH
//
//  Created by Charlie Close on 17/03/2025.
//

#ifndef eos_h
#define eos_h

#include <metal_stdlib>
using namespace metal;

static inline float4 eos(int materialId,
           float u,
           float density,
           texture2d<float, access::sample> texIron,
           texture2d<float, access::sample> texForesite,
           texture2d<float, access::sample> texFe,
           texture2d<float, access::sample> texHHe,
           texture2d<float, access::sample> texIce,
           texture2d<float, access::sample> texRock
           ) {
    texture2d<float, access::sample> texture;
    float2 uv;
    sampler textureSampler(mag_filter::linear, min_filter::linear, mip_filter::none, address::clamp_to_edge);
    switch (materialId) {
        case (200): {
            float uCoord = log(u / 1e-8) / log(5e9 / 1e4);
            float rhoCoord = log(density / 1e-7 ) / log(1000 / 0.1);
            uv = clamp(float2(uCoord, rhoCoord), float2(0.0), float2(1.0));
            texture = texHHe;
            break;
        }
        case (201): {
            float uCoord = log(u / 1e-9) / log(5e9 / 1e3);
            float rhoCoord = log(density / 1e-6 ) / log(6000 / 1.);
            uv = clamp(float2(uCoord, rhoCoord), float2(0.0), float2(1.0));
            texture = texIce;
            break;
        }
        case (202): {
            float uCoord = log(u / 1e-8) / log(1e9 / 1e4);
            float rhoCoord = log(density / 1e-6 ) / log(20000 / 1.);
            uv = clamp(float2(uCoord, rhoCoord), float2(0.0), float2(1.0));
            texture = texRock;
            break;
        }
        case (400): {
            float uCoord = log(u / ANEOS_MIN_U) / log(ANEOS_MAX_U / ANEOS_MIN_U);
            float rhoCoord = log(density / ANEOS_MIN_RHO) / log(ANEOS_MAX_RHO / ANEOS_MIN_RHO);
            uv = clamp(float2(uCoord, rhoCoord), float2(0.0), float2(1.0));
            texture = texForesite;
            break;
        }
        case (401): {
            float uCoord = log(u / ANEOS_MIN_U) / log(ANEOS_MAX_U / ANEOS_MIN_U);
            float rhoCoord = log(density / ANEOS_MIN_RHO) / log(ANEOS_MAX_RHO / ANEOS_MIN_RHO);
            uv = clamp(float2(uCoord, rhoCoord), float2(0.0), float2(1.0));
            texture = texIron;
            break;
        }
        case(402): {
            float uCoord = log(u / ANEOS_MIN_U) / log(ANEOS_MAX_U / ANEOS_MIN_U);
            float rhoCoord = log(density / ANEOS_MIN_RHO) / log(ANEOS_MAX_RHO / ANEOS_MIN_RHO);
            uv = clamp(float2(uCoord, rhoCoord), float2(0.0), float2(1.0));
            texture = texFe;
            break;
        }
    }
    return texture.sample(textureSampler, uv);
}


#endif /* eos_h */
