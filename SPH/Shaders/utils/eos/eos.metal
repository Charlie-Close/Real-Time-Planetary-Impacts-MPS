//
//  eos.metal
//  SPH
//
//  Created by Charlie Close on 03/11/2024.
//

#include "eos.h"
#include "materials.h"

float calculateP_c(TillotsonEOS material, float eta, float w_inv, float density, float internalEnergy, float mu) {
    float P_c = 0;
    if (eta < material.eta_zero) {
        P_c = 0;
    } else {
        P_c = (material.a + material.b * w_inv) * density * internalEnergy + material.A * mu + material.B * mu * mu;
        if (eta < material.eta_min) {
            P_c *= (eta - material.eta_zero) / (material.eta_min - material.eta_zero);
        }
    }
    return P_c;
}

float calculateP_e(TillotsonEOS material, float density, float internalEnergy, float w_inv, float mu, float nu) {
    return material.a * density * internalEnergy + (material.b * density * internalEnergy * w_inv + material.A * mu * fast::exp( - material.beta * nu )) * fast::exp( - material.alpha * nu * nu );
}

TillotsonEOS materialFromID(int id) {
    switch (id) {
        case 100:
            return TillotsonIron;
            break;
        case 101:
            return TillotsonGranite;
            break;
        case 102:
            return TillotsonWater;
            break;
        case 103:
            return TillotsonBasalt;
            break;
            
            
        case 402:
            return TillotsonIron;
            break;
        case 400:
            return TillotsonGranite;
            break;
        case 404:
            return TillotsonWater;
            break;
        case 401:
            return TillotsonBasalt;
            break;
        default:
            return TillotsonBasalt;
            break;
    }
}

float pressureEos(float density, float internalEnergy, int materialId) {
    TillotsonEOS material = materialFromID(materialId);
    
    float eta = density / material.rho_0;
    float eta_sq = eta * eta;
    float mu = eta - 1;
    float nu = 1 / eta - 1;
    float w = internalEnergy / (material.u_0 * eta_sq) + 1;
    float w_inv = 1 / w;
    
    float P = 0;
    
    if (1.0 < eta || internalEnergy < material.u_iv) {
        // Cold state
        P = calculateP_c(material, eta, w_inv, density, internalEnergy, mu);
    } else if (eta < 1 && material.u_cv < internalEnergy) {
        // Hot state
        P = calculateP_e(material, density, internalEnergy, w_inv, mu, nu);
    } else {
        // Hybrid state
        float P_c = calculateP_c(material, eta, w_inv, density, internalEnergy, mu);
        float P_e = calculateP_e(material, density, internalEnergy, w_inv, mu, nu);
        P = ((internalEnergy - material.u_iv) * P_e + (material.u_cv - internalEnergy) * P_c) / (material.u_cv - material.u_iv);
    }
    
    if (P < material.P_min) {
        P = material.P_min;
    }
    
    return P;
}

float speedOfSoundEos(float density, float internalEnergy, int materialId) {
    TillotsonEOS material = materialFromID(materialId);
    
    const float rho_0_inv = 1.0f / material.rho_0;
    const float eta = density * rho_0_inv;
    const float rho_inv = 1.0f / density;
    const float eta_sq = eta * eta;
    const float mu = eta - 1.0f;
    const float nu = 1.0f / eta - 1.0f;
    const float w = internalEnergy / (material.u_0 * eta_sq) + 1.0f;
    const float w_inv = 1.0f / w;
    const float w_inv_sq = w_inv * w_inv;
    const float exp_beta = fast::exp(-material.beta * nu);
    const float exp_alpha = fast::exp(-material.alpha * nu * nu);

    // Pressure and sound speed squared for the cold state
    float P_c = (material.a + material.b * w_inv) * density * internalEnergy + material.A * mu + material.B * mu * mu;
    if (eta < material.eta_zero) {
        P_c = 0.0f;
    } else if (eta < material.eta_min) {
        P_c *= (eta - material.eta_zero) / (material.eta_min - material.eta_zero);
    }
    float c_sq_c = P_c * rho_inv * (1.0f + material.a + material.b * w_inv)
                   + material.b * (w - 1.0f) * w_inv_sq * (2.0f * internalEnergy - P_c * rho_inv)
                   + rho_inv * (material.A + material.B * (eta_sq - 1.0f));

    // Pressure and sound speed squared for the expanded state
    float P_e = material.a * density * internalEnergy + (material.b * density * internalEnergy * w_inv + material.A * mu * exp_beta) * exp_alpha;
    float c_sq_e = P_e * rho_inv * (1.0f + material.a + material.b * w_inv * exp_alpha)
                   + (material.b * density * internalEnergy * w_inv_sq / eta_sq *
                      (rho_inv / material.u_0 * (2.0f * internalEnergy - P_e * rho_inv)
                       + 2.0f * material.alpha * nu * w * rho_0_inv)
                      + material.A * rho_0_inv * (1.0f + mu / eta_sq * (material.beta + 2.0f * material.alpha * nu - eta))
                      * exp_beta) * exp_alpha;

    // Choose the correct sound speed squared value based on the state
    float c_sq;
    if (eta > 1.0f || internalEnergy < material.u_iv) {
        // Cold state
        c_sq = c_sq_c;
    } else if (eta < 1.0f && internalEnergy > material.u_cv) {
        // Hot state
        c_sq = c_sq_e;
    } else {
        // Hybrid state (interpolation between cold and hot)
        c_sq = ((internalEnergy - material.u_iv) * c_sq_e + (material.u_cv - internalEnergy) * c_sq_c) / (material.u_cv - material.u_iv);
    }

    // Ensure sound speed squared is not less than a minimum threshold
    c_sq = fmax(c_sq, material.A * rho_0_inv);

    // Return the sound speed as the square root of c_sq
    return fast::sqrt(c_sq);
}
