//
//  kernels.h
//  SPH
//
//  Created by Charlie Close on 13/10/2024.
//

#ifndef kernels_h
#define kernels_h

#include <metal_stdlib>
#include "../../Parameters.h"
using namespace metal;



constant float gamma1 = 1 / GAMMA;

static inline float W(float r, float h1) {
#ifdef WENDLAND_C2_KERNEL
    float H1 = h1 * gamma1;
    float u = r * H1;
    if (u > 1.f) {
        return 0.f;
    }
    float u_2 = u * u;
    float u_3 = u_2 * u;
    float u_4 = u_3 * u;
    float u_5 = u_4 * u;
    float fac = KERNEL_CONSTANT * H1 * H1 * H1;

    return fac * (4 * u_5 - 15 * u_4 + 20 * u_3 - 10 * u_2 + 1);
#elifdef CUBIC_SPLINE_KERNEL
    float H1 = h1 * gamma1;
    float u = r * H1;
    if (u > 1.f) {
        return 0.f;
    }
    
    float u_2 = u * u;
    float u_3 = u_2 * u;
    float fac = KERNEL_CONSTANT * H1 * H1 * H1;
    
    if (u > 0.5) {
//        −u3 + 3u2 − 3u + 1
        return fac * (- u_3 + 3 * u_2 - 3 * u + 1);
    }
//    3u3 − 3u2 + 1/2
    return fac * (3 * u_3 - 3 * u_2 + 0.5);
#elifdef QUARTIC_SPLINE_KERNEL
    float H1 = h1 * gamma1;
    float u = r * H1;
    if (u > 1.f) {
        return 0.f;
    }
    
    float u_2 = u * u;
    float u_3 = u_2 * u;
    float u_4 = u_3 * u;
    float fac = KERNEL_CONSTANT * H1 * H1 * H1;
    
    if (u > 0.6) {
//        u4 − 4u3 + 6u2 − 4u + 1
        return fac * (u_4 - 4 * u_3 + 6 * u_2 - 4 * u + 1);
    } else if (u > 0.2) {
//        −4u4 + 8u3 − 24/5 u2 + 8/25 u + 44/125
        return fac * (- 4 * u_4 + 8 * u_3 - 4.8 * u_2 + 0.32 * u + 0.352);
    }
// 6u4 − 12/5 u2 + 46/125
    return fac * (6 * u_4 - 2.4 * u_2 + 0.368);


#endif
}

static inline float3 gradW(float3 x_ij, float r, float r1, float h1) {
#ifdef WENDLAND_C2_KERNEL
    float H1 = h1 * gamma1;
    float u = r * H1;
    if (u > 1) {
        return { 0, 0, 0 };
    }
    float u_2 = u * u;
    float u_3 = u_2 * u;
    float u_4 = u_3 * u;
    float H1_2 = H1 * H1;
    float H1_4 = H1_2 * H1_2;
    float fac = KERNEL_CONSTANT * H1_4;
    float3 r_hat = x_ij * r1;
    
    return fac * (20 * u_4 - 60 * u_3 + 60 * u_2 - 20 * u) * r_hat;
#elifdef CUBIC_SPLINE_KERNEL
    float H1 = h1 * gamma1;
    float u = r * H1;
    if (u > 1) {
        return { 0, 0, 0 };
    }
    float u_2 = u * u;
    float H1_2 = H1 * H1;
    float H1_4 = H1_2 * H1_2;
    float fac = KERNEL_CONSTANT * H1_4;
    float3 r_hat = x_ij * r1;
    
    if (u > 0.5) {
//        −u3 + 3u2 − 3u + 1
        return fac * (- 3 * u_2 + 6 * u - 3) * r_hat;
    }
//    3u3 − 3u2 + 1/2
    return fac * (9 * u_2 - 6 * u) * r_hat;
#elifdef QUARTIC_SPLINE_KERNEL
    float H1 = h1 * gamma1;
    float u = r * H1;
    if (u > 1) {
        return { 0, 0, 0 };
    }
    float u_2 = u * u;
    float u_3 = u_2 * u;
    float H1_2 = H1 * H1;
    float H1_4 = H1_2 * H1_2;
    float fac = KERNEL_CONSTANT * H1_4;
    float3 r_hat = x_ij * r1;
    
    if (u > 0.6) {
//        u4 − 4u3 + 6u2 − 4u + 1
        return fac * (4 * u_3 - 12 * u_2 + 12 * u - 4) * r_hat;
    } else if (u > 0.2) {
//        −4u4 + 8u3 − 24/5 u2 + 8/25 u + 44/125
        return fac * (-16 * u_3 + 24 * u_2 - 9.6 * u + 0.32) * r_hat;
    }
// 6u4 − 12/5 u2 + 46/125
    return fac * (24 * u_3 - 4.8 * u) * r_hat;


#endif
}

static inline float dW_dh(float r, float h1) {
#ifdef WENDLAND_C2_KERNEL
    float H1 = h1 * gamma1;
    float u = r * H1;
    if (u > 1) {
        return 0;
    }
    float u_2 = u * u;
    float u_3 = u_2 * u;
    float u_4 = u_3 * u;
    float u_5 = u_4 * u;
    float H1_2 = H1 * H1;
    float H1_4 = H1_2 * H1_2;
    float fac = KERNEL_CONSTANT * H1_4;

    return - fac * (32 * u_5 - 105 * u_4 + 120 * u_3 - 50 * u_2 + 3);
#elifdef CUBIC_SPLINE_KERNEL
    float H1 = h1 * gamma1;
    float u = r * H1;
    if (u > 1) {
        return 0;
    }
    float u_2 = u * u;
    float u_3 = u_2 * u;
    float H1_2 = H1 * H1;
    float H1_4 = H1_2 * H1_2;
    float fac = KERNEL_CONSTANT * H1_4;

    if (u > 0.5) {
//        3f = −3u3 + 9u2 − 9u + 3
//        uf' =  - 3 * u_3 + 6 * u_2 - 3u;
        return - fac * ( - 6 * u_3 + 15 * u_2 - 12 * u + 3);
    }
//    3f = 9u3 − 9u2 + 1.5
//    uf' = 9 * u_3 - 6 * u_2
    return - fac * ( 18 * u_3 - 15 * u_2 + 1.5);
#elifdef QUARTIC_SPLINE_KERNEL
    float H1 = h1 * gamma1;
    float u = r * H1;
    if (u > 1) {
        return 0;
    }
    float u_2 = u * u;
    float u_3 = u_2 * u;
    float u_4 = u_3 * u;
    float H1_2 = H1 * H1;
    float H1_4 = H1_2 * H1_2;
    float fac = KERNEL_CONSTANT * H1_4;
    
    if (u > 0.6) {
//        3f = 3u4 − 12u3 + 18u2 − 12u + 3
//      uf'  = 4u4 - 12u3 + 12u2 - 4u
        return fac * ( 7 * u_4 - 24 * u_3 + 30 * u_2 - 16 * u + 3);
//        return fac * (4 * u_3 - 12 * u_2 + 12 * u - 4) * r_hat;
    } else if (u > 0.2) {
//        3f = −12u4 + 24u3 − 14.4u2 + 0.96u + 1.056
//       uf' = -16u4 + 24u3 - 9.6u2 + 0.32u
        return fac * (-28 * u_4 + 48 * u_3 - 24 * u_2 + 1.28 * u + 1.056);
    }
// 3f = 18u4 − 7.2u2 + 1.104
// uf'= 24u4 - 4.8u2
    return fac * (32 * u_4 - 12 * u_2 + 1.104);
#endif
}

static inline float dphi_dr(float3 x_ij, float eta) {
    float eta1 = 1.f / eta;
    float r = fast::length(x_ij);
    float q = r * eta1;
    
    if (q <= 1) {
        return (eta1 * eta1) * ((4 / 3) * q - (6 / 5) * q * q * q + 0.5 * q * q * q * q);
    }
    if (q <= 2) {
        float q1 = 1.f / q;
        return (eta1 * eta1) * ((8 / 3) * q - 3 * q * q + (6 / 5) * q * q * q - (1 / 6) * q * q * q * q - (1 / 15) * q1 * q1);
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
