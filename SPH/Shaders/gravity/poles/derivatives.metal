//
//  derivatives.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include "poles.h"
#include <metal_stdlib>
using namespace metal;

float D_soft_1(const float u) {
    /* -3u^7 + 15u^6 - 28u^5 + 21u^4 - 7u^2 + 3 */
    float phi = -3.f * u + 15.f;
    phi = phi * u - 28.f;
    phi = phi * u + 21.f;
    phi = phi * u;
    phi = phi * u - 7.f;
    phi = phi * u;
    phi = phi * u + 3.f;

    return phi;
}

float D_soft_2(const float u) {
    /* -21u^6 + 90u^5 - 140u^4 + 84u^3 - 14u */
    float phi = -21.f * u + 90.f;
    phi = phi * u - 140.f;
    phi = phi * u + 84.f;
    phi = phi * u;
    phi = phi * u - 14.f;
    phi = phi * u;

    return phi;
}

float D_soft_3(const float u) {
    /* -105u^5 + 360u^4 - 420u^3 + 168u^2 */
    float phi = -105.f * u + 360.f;
    phi = phi * u - 420.f;
    phi = phi * u + 168.f;
    phi = phi * u;
    phi = phi * u;

    return phi;
}

float D_soft_4(const float u) {
    /* -315u^4 + 720u^3 - 420u^2 */
    float phi = -315.f * u + 720.f;
    phi = phi * u - 420.f;
    phi = phi * u;
    phi = phi * u;

    return phi;
}

float D_soft_5(const float u) {
  /* -315u^3 + 420u */
  float phi = -315.f * u;
  phi = phi * u + 420.f;
  phi = phi * u;

  return phi;
}

Derivatives derivatives(float3 vec, float eps) {
    float Dt_1;
#if P > 0
    float Dt_2;
#endif
#if P > 1
    float Dt_3;
#endif
#if P > 2
    float Dt_4;
#endif
#if P > 3
    float Dt_5;
#endif
    float r2 = length_squared(vec);
    float r = sqrt(r2);
    float r_inv = r > 1e-6 ? 1 / r : 0;
    float eps2 = eps * eps;
    
    if (r2 < eps2) {
        const float eps_inv = 1.f / eps;
        const float u = r * eps_inv;

        Dt_1 = eps_inv * D_soft_1(u);
#if P > 0
        const float eps_inv2 = eps_inv * eps_inv;
        Dt_2 = eps_inv2 * D_soft_2(u);
#endif
#if P > 1
        const float eps_inv3 = eps_inv2 * eps_inv;
        Dt_3 = eps_inv3 * D_soft_3(u);
#endif
#if P > 2
        const float eps_inv4 = eps_inv3 * eps_inv;
        Dt_4 = eps_inv4 * D_soft_4(u);
#endif
#if P > 3
        const float eps_inv5 = eps_inv4 * eps_inv;
        Dt_5 = eps_inv5 * D_soft_5(u);
#endif
    } else {
        Dt_1 = r_inv; /* 1 / r */
#if P > 0
        Dt_2 = -1.f * Dt_1 * r_inv; /* -1 / r^2 */
#endif
#if P > 1
        Dt_3 = -3.f * Dt_2 * r_inv; /* 3 / r^3 */
#endif
#if P > 2
        Dt_4 = -5.f * Dt_3 * r_inv; /* -15 / r^4 */
#endif
#if P > 3
        Dt_5 = -7.f * Dt_4 * r_inv; /* 105 / r^5 */
#endif
    }
    
#if P > 0
    const float3 r_r = vec * r_inv;
#endif
#if P > 1
    const float3 r_r2 = r_r * r_r;
#endif
#if P > 2
    const float3 r_r3 = r_r2 * r_r;
#endif
#if P > 3
    const float3 r_r4 = r_r3 * r_r;
#endif
    
    Derivatives dev;
    
    dev.expansion[M] = Dt_1;
#if P > 0
    dev.expansion[X] = r_r.x * Dt_2;
    dev.expansion[Y] = r_r.y * Dt_2;
    dev.expansion[Z] = r_r.z * Dt_2;
#endif
#if P > 1
    Dt_2 *= r_inv;

    dev.expansion[XX] = r_r2.x * Dt_3 + Dt_2;
    dev.expansion[XY] = r_r.x * r_r.y * Dt_3;
    dev.expansion[XZ] = r_r.x * r_r.z * Dt_3;
    dev.expansion[YY] = r_r2.y * Dt_3 + Dt_2;
    dev.expansion[YZ] = r_r.y * r_r.z * Dt_3;
    dev.expansion[ZZ] = r_r2.z * Dt_3 + Dt_2;
    
#endif
#if P > 2
    Dt_3 *= r_inv;

    /* 3rd order derivatives */
    dev.expansion[XXX] = r_r3.x * Dt_4 + 3.f * r_r.x * Dt_3;
    dev.expansion[YYY] = r_r3.y * Dt_4 + 3.f * r_r.y * Dt_3;
    dev.expansion[ZZZ] = r_r3.z * Dt_4 + 3.f * r_r.z * Dt_3;
    dev.expansion[XXY] = r_r2.x * r_r.y * Dt_4 + r_r.y * Dt_3;
    dev.expansion[XXZ] = r_r2.x * r_r.z * Dt_4 + r_r.z * Dt_3;
    dev.expansion[XYY] = r_r2.y * r_r.x * Dt_4 + r_r.x * Dt_3;
    dev.expansion[YYZ] = r_r2.y * r_r.z * Dt_4 + r_r.z * Dt_3;
    dev.expansion[XZZ] = r_r2.z * r_r.x * Dt_4 + r_r.x * Dt_3;
    dev.expansion[YZZ] = r_r2.z * r_r.y * Dt_4 + r_r.y * Dt_3;
    dev.expansion[XYZ] = r_r.x * r_r.y * r_r.z * Dt_4;
#endif
#if P > 3
      Dt_3 *= r_inv;
      Dt_4 *= r_inv;

      /* 4th order derivatives */
      dev.expansion[XXXX] = r_r4.x * Dt_5 + 6.f * r_r2.x * Dt_4 + 3.f * Dt_3;
      dev.expansion[YYYY] = r_r4.y * Dt_5 + 6.f * r_r2.y * Dt_4 + 3.f * Dt_3;
      dev.expansion[ZZZZ] = r_r4.z * Dt_5 + 6.f * r_r2.z * Dt_4 + 3.f * Dt_3;
      dev.expansion[XXXY] = r_r3.x * r_r.y * Dt_5 + 3.f * r_r.x * r_r.y * Dt_4;
      dev.expansion[XXXZ] = r_r3.x * r_r.z * Dt_5 + 3.f * r_r.x * r_r.z * Dt_4;
      dev.expansion[XYYY] = r_r3.y * r_r.x * Dt_5 + 3.f * r_r.y * r_r.x * Dt_4;
      dev.expansion[YYYZ] = r_r3.y * r_r.z * Dt_5 + 3.f * r_r.y * r_r.z * Dt_4;
      dev.expansion[XZZZ] = r_r3.z * r_r.x * Dt_5 + 3.f * r_r.z * r_r.x * Dt_4;
      dev.expansion[YZZZ] = r_r3.z * r_r.y * Dt_5 + 3.f * r_r.z * r_r.y * Dt_4;
      dev.expansion[XXYY] = r_r2.x * r_r2.y * Dt_5 + r_r2.x * Dt_4 + r_r2.y * Dt_4 + Dt_3;
      dev.expansion[XXZZ] = r_r2.x * r_r2.z * Dt_5 + r_r2.x * Dt_4 + r_r2.z * Dt_4 + Dt_3;
      dev.expansion[YYZZ] = r_r2.y * r_r2.z * Dt_5 + r_r2.y * Dt_4 + r_r2.z * Dt_4 + Dt_3;
      dev.expansion[XXYZ] = r_r2.x * r_r.y * r_r.z * Dt_5 + r_r.y * r_r.z * Dt_4;
      dev.expansion[XYYZ] = r_r2.y * r_r.x * r_r.z * Dt_5 + r_r.x * r_r.z * Dt_4;
      dev.expansion[XYZZ] = r_r2.z * r_r.x * r_r.y * Dt_5 + r_r.x * r_r.y * Dt_4;
#endif
    return dev;
}
