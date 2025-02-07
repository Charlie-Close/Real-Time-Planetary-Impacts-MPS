//
//  accept.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include "poles.h"
#include <metal_stdlib>
#include "../../../Parameters.h"
using namespace metal;


bool gravity_M2L_accept(Multipole A, Multipole B) {
    /* Order of the expansion */
    const int p = P;

    /* Sizes of the multipoles */
    const float rho_A = A.size;
    const float rho_B = B.size;
    const float r2 = length_squared(A.pos - B.pos);

    /* Max size of both multipoles */
    const float rho_max = max(rho_A, rho_B);

    /* Compute the error estimator (without the 1/M_B term that cancels out) */
    float E_BA_term = 0.f;
    for (int n = 0; n <= p; ++n) {
    E_BA_term += binomial(p, n) * B.power[n] * integer_powf(rho_A, p - n);
    }
    E_BA_term *= 8.f;
    if (rho_A + rho_B > 0.f) {
        E_BA_term *= rho_max;
        E_BA_term /= (rho_A + rho_B);
    }
    
    /* Compute r^p = (r^2)^(p/2) */
    const float r_to_p = integer_powf(r2, (p / 2));

    float f_MAC_inv = r2;

    /* Get the mimimal acceleration in A */
    const float min_a_grav = A.minGrav;

    /* Get the sum of the multipole sizes */
    const float rho_sum = rho_A + rho_B;

    /* Test the different conditions */
    /* Condition 1: We are in the converging part of the Taylor expansion */
    const int cond_1 = rho_sum * rho_sum < r2;

    /* Condition 2: The contribution is accurate enough
     * (E_BA * (1 / r^(p)) * ((1 / r^2) * W) < eps * a_min) */
    const int cond_2 = E_BA_term < GRAVITY_ETA * min_a_grav * r_to_p * f_MAC_inv;

    return cond_1 && cond_2;
}

