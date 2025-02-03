//
//  materials.h
//  SPH
//
//  Created by Charlie Close on 03/11/2024.
//

#ifndef materials_h
#define materials_h

struct TillotsonEOS {
    float rho_0;     // Reference density
    float a;         // Dimensionless coefficient for pressure
    float b;         // Dimensionless coefficient for pressure
    float A;         // Coefficient for condensed phase pressure
    float B;         // Coefficient for condensed phase pressure
    float u_0;       // Specific internal energy where vaporization starts
    float u_iv;      // Specific internal energy for complete vaporization
    float u_cv;      // Specific internal energy for complete vaporization (alternative term)
    float alpha;     // Coefficient for pressure in vapor phase
    float beta;      // Coefficient affecting high-energy behavior
    float eta_min;   // Minimum eta value for this material
    float P_min;     // Minimum pressure
    float eta_zero;  // Eta zero value for this material
};

// Iron
constant TillotsonEOS TillotsonIron = {
    7800.0f,          // rho_0 (7800 / 1000)
    0.5f,          // a
    1.5f,          // b
    1.28e11f,      // A (1.28e11 / 1e16)
    1.05e11f,      // B (1.05e11 / 1e16)
    9.5e6f,      // u_0 (9.5e6 / 1e16)
    2.4e6f,      // u_iv (2.4e6 / 1e16)
    8.67e6f,     // u_cv (8.67e6 / 1e16)
    5.0f,          // alpha
    5.0f,          // beta
    0.0f,          // eta_min
    0.0f,          // P_min
    0.0f           // eta_zero
};

// Granite
constant TillotsonEOS TillotsonGranite = {
    2680.f,         // rho_0 (2680 / 1000)
    0.5f,          // a
    1.3f,          // b
    1.8e10f,      // A (1.8e10 / 1e16)
    1.8e10f,      // B (1.8e10 / 1e16)
    1.6e7f,      // u_0 (1.6e7 / 1e16)
    3.5e6f,     // u_iv (3.5e6 / 1e16)
    1.8e7f,      // u_cv (1.8e7 / 1e16)
    5.0f,          // alpha
    5.0f,          // beta
    0.0f,          // eta_min
    0.0f,          // P_min
    0.0f           // eta_zero
};

// Basalt
constant TillotsonEOS TillotsonBasalt = {
    2700.f,          // rho_0 (2700 / 1000)
    0.5f,          // a
    1.5f,          // b
    2.67e10f,     // A (2.67e10 / 1e16)
    2.67e10f,     // B (2.67e10 / 1e16)
    4.87e8f,    // u_0 (4.87e8 / 1e16)
    4.72e6f,   // u_iv (4.72e6 / 1e16)
    1.82e7f,    // u_cv (1.82e7 / 1e16)
    5.0f,          // alpha
    5.0f,          // beta
    0.0f,          // eta_min
    0.0f,          // P_min
    0.0f           // eta_zero
};

// Water
constant TillotsonEOS TillotsonWater = {
    998.f,        // rho_0 (998 / 1000)
    0.7f,          // a
    0.15f,         // b
    2.18e9f,    // A (2.18e9 / 1e16)
    1.325e10f,    // B (1.325e10 / 1e16)
    7.0e6f,     // u_0 (7.0e6 / 1e16)
    4.19e5f,  // u_iv (4.19e5 / 1e16)
    2.69e6f,   // u_cv (2.69e6 / 1e16)
    10.0f,         // alpha
    5.0f,          // beta
    0.925f,        // eta_min
    0.0f,          // P_min
    0.875f         // eta_zero
};


//// Iron
//constant TillotsonEOS TillotsonIron = {
//    7.8f,          // rho_0 (7800 / 1000)
//    0.5f,          // a
//    1.5f,          // b
//    1.28e-5f,      // A (1.28e11 / 1e16)
//    1.05e-5f,      // B (1.05e11 / 1e16)
//    0.95e-9f,      // u_0 (9.5e6 / 1e16)
//    0.24e-9f,      // u_iv (2.4e6 / 1e16)
//    0.867e-9f,     // u_cv (8.67e6 / 1e16)
//    5.0f,          // alpha
//    5.0f,          // beta
//    0.0f,          // eta_min
//    0.0f,          // P_min
//    0.0f           // eta_zero
//};
//
//// Granite
//constant TillotsonEOS TillotsonGranite = {
//    2.68f,         // rho_0 (2680 / 1000)
//    0.5f,          // a
//    1.3f,          // b
//    0.18e-5f,      // A (1.8e10 / 1e16)
//    0.18e-5f,      // B (1.8e10 / 1e16)
//    0.16e-9f,      // u_0 (1.6e7 / 1e16)
//    0.035e-9f,     // u_iv (3.5e6 / 1e16)
//    0.18e-9f,      // u_cv (1.8e7 / 1e16)
//    5.0f,          // alpha
//    5.0f,          // beta
//    0.0f,          // eta_min
//    0.0f,          // P_min
//    0.0f           // eta_zero
//};
//
//// Basalt
//constant TillotsonEOS TillotsonBasalt = {
//    2.7f,          // rho_0 (2700 / 1000)
//    0.5f,          // a
//    1.5f,          // b
//    0.267e-5f,     // A (2.67e10 / 1e16)
//    0.267e-5f,     // B (2.67e10 / 1e16)
//    0.0487e-9f,    // u_0 (4.87e8 / 1e16)
//    0.00472e-9f,   // u_iv (4.72e6 / 1e16)
//    0.0182e-9f,    // u_cv (1.82e7 / 1e16)
//    5.0f,          // alpha
//    5.0f,          // beta
//    0.0f,          // eta_min
//    0.0f,          // P_min
//    0.0f           // eta_zero
//};
//
//// Water
//constant TillotsonEOS TillotsonWater = {
//    0.998f,        // rho_0 (998 / 1000)
//    0.7f,          // a
//    0.15f,         // b
//    0.0218e-5f,    // A (2.18e9 / 1e16)
//    0.1325e-5f,    // B (1.325e10 / 1e16)
//    0.007e-9f,     // u_0 (7.0e6 / 1e16)
//    0.000419e-9f,  // u_iv (4.19e5 / 1e16)
//    0.00269e-9f,   // u_cv (2.69e6 / 1e16)
//    10.0f,         // alpha
//    5.0f,          // beta
//    0.925f,        // eta_min
//    0.0f,          // P_min
//    0.875f         // eta_zero
//};

#endif /* materials_h */
