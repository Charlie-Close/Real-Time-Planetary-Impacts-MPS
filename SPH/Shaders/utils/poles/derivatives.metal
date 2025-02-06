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
    }
    
#if P > 0
    const float rx_r = vec.x * r_inv;
    const float ry_r = vec.y * r_inv;
    const float rz_r = vec.z * r_inv;
#endif
#if P > 1
    const float rx_r2 = rx_r * rx_r;
    const float ry_r2 = ry_r * ry_r;
    const float rz_r2 = rz_r * rz_r;
#endif
#if P > 2
    const float rx_r3 = rx_r2 * rx_r;
    const float ry_r3 = ry_r2 * ry_r;
    const float rz_r3 = rz_r2 * rz_r;
#endif
    
    Derivatives dev;
    
    dev.expansion[M] = Dt_1;
#if P > 0
    dev.expansion[X] = rx_r * Dt_2;
    dev.expansion[Y] = ry_r * Dt_2;
    dev.expansion[Z] = rz_r * Dt_2;
#endif
#if P > 1
    Dt_2 *= r_inv;
    
    dev.expansion[XX] = rx_r2 * Dt_3 + Dt_2;
    dev.expansion[XY] = rx_r * ry_r * Dt_3;
    dev.expansion[XZ] = rx_r * rz_r * Dt_3;
    dev.expansion[YY] = ry_r2 * Dt_3 + Dt_2;
    dev.expansion[YZ] = ry_r * rz_r * Dt_3;
    dev.expansion[ZZ] = rz_r2 * Dt_3 + Dt_2;
    
#endif
#if P > 2
    Dt_3 *= r_inv;

    /* 3rd order derivatives */
    dev.expansion[XXX] = rx_r3 * Dt_4 + 3.f * rx_r * Dt_3;
    dev.expansion[YYY] = ry_r3 * Dt_4 + 3.f * ry_r * Dt_3;
    dev.expansion[ZZZ] = rz_r3 * Dt_4 + 3.f * rz_r * Dt_3;
    dev.expansion[XXY] = rx_r2 * ry_r * Dt_4 + ry_r * Dt_3;
    dev.expansion[XXZ] = rx_r2 * rz_r * Dt_4 + rz_r * Dt_3;
    dev.expansion[XYY] = ry_r2 * rx_r * Dt_4 + rx_r * Dt_3;
    dev.expansion[YYZ] = ry_r2 * rz_r * Dt_4 + rz_r * Dt_3;
    dev.expansion[XZZ] = rz_r2 * rx_r * Dt_4 + rx_r * Dt_3;
    dev.expansion[YZZ] = rz_r2 * ry_r * Dt_4 + ry_r * Dt_3;
    dev.expansion[XYZ] = rx_r * ry_r * rz_r * Dt_4;
#endif
    return dev;
}
