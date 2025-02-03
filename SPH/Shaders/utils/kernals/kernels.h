//
//  kernels.h
//  SPH
//
//  Created by Charlie Close on 13/10/2024.
//

#ifndef kernels_h
#define kernels_h

#include <metal_stdlib>
using namespace metal;

constant float pi1 = 1 / M_PI_F;

static inline float W(float3 x_ij, float h) {
    float r = fast::length(x_ij);
    float h1 = 1 / h;
    float q = r * h1;
    
    if (q > 2.0) {
        return 0;
    }
    
    float fac = pi1 * h1 * h1 * h1;
    
    if (q > 1.0) {
        float tmp = 2.f - q;
        return fac * 0.25 * tmp * tmp * tmp;
    }
    
    return fac * (1 - 1.5 * q * q * (1 - 0.5 * q));
}

static inline float3 gradW(float3 x_ij, float h) {
    float h1 = 1.f / h;
    float r = fast::length(x_ij);
    float q = r * h1;
    
    if (q > 2.f or q < 1e-8) {
        return { 0, 0, 0 };
    }
    
    float3 r_hat = x_ij / r;
    
    float fac = pi1 * h1 * h1 * h1 * h1;
    
    if (q > 1.f) {
        float tmp = 2.f - q;
        return -0.75 * fac * tmp * tmp * r_hat;
    }
    
    return 3 * fac * (0.75 * q - 1) * q  * r_hat;
}

static inline float dW_dh(float3 x_ij, float h) {
    float h1 = 1.f / h;
    float r = fast::length(x_ij);
    float q = r * h1;
    
    if (q > 2.0) {
        return 0;
    }
    
    float fac = - 3 * pi1 * h1 * h1 * h1 * h1;
    
    if (q > 1.0) {
        return fac * (2 - 4 * q + 2.5 * q * q - 0.5 * q * q * q);
    }
    
    return fac * (1 + 0.5 * (3 * q - 5) * q * q);
}

static inline float dphi_dr(float3 x_ij, float h) {
    float h1 = 1.f / min(h, 0.16);
    float r = fast::length(x_ij);
    float q = r * h1;
    
    if (q <= 1) {
        return (h1 * h1) * ((4 / 3) * q - (6 / 5) * q * q * q + 0.5 * q * q * q * q);
    }
    if (q <= 2) {
        float q1 = 1.f / q;
        return (h1 * h1) * ((8 / 3) * q - 3 * q * q + (6 / 5) * q * q * q - (1 / 6) * q * q * q * q - (1 / 15) * q1 * q1);
    }
    
    float r1 = 1.f / r;
    return r1 * r1;
}
static inline float dphi_dh(float3 x_ij, float h) {
    float h1 = 1.f / min(h, 0.16);
    float r = fast::length(x_ij);
    float q = r * h1;
    
    if (q >= 2) {
        return 0;
    }
    
    float fact = h1 * h1;
    
    if (r > 1) {
        return fact * (-2 * q * q + (6 / 5) * q * q * q - (3 / 5) * q * q * q * q * q + (7 / 5));
    }
    
    return fact * (-4 * q * q + 4 * q * q * q - 1.5 * q * q * q * q + (1 / 5) * q * q * q * q * q + (8 / 5));
}

#endif /* kernels_h */
