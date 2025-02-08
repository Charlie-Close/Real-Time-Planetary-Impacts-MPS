//
//  poles.h
//  SPH
//
//  Created by Charlie Close on 23/01/2025.
//

#ifndef poles_h
#define poles_h


#include <metal_stdlib>
using namespace metal;
#include "../../../Parameters.h"

constant float G = 6.67e-5;

constant uint M = 0;
constant uint X = 1;
constant uint Y = 2;
constant uint Z = 3;

constant uint XX = 4;
constant uint XY = 5;
constant uint XZ = 6;
constant uint YY = 7;
constant uint YZ = 8;
constant uint ZZ = 9;

constant uint XXX = 10;
constant uint XXY = 11;
constant uint XXZ = 12;
constant uint XYY = 13;
constant uint XYZ = 14;
constant uint XZZ = 15;
constant uint YYY = 16;
constant uint YYZ = 17;
constant uint YZZ = 18;
constant uint ZZZ = 19;


// Q_m
// Q_x, Q_y, Q_z
// Q_xx, Q_xy, Q_xz, Q_yy, Q_yz, Q_zz

typedef struct {
    float3 pos;
    float3 min;
    float3 max;
    float size;
    float expansion[N_EXPANSION_TERMS];
    float power[P+1];
    float minGrav;
    bool active;
} Multipole;


typedef struct {
    float3 pos;
    float expansion[N_EXPANSION_TERMS];
} Local;

typedef struct {
    float expansion[N_EXPANSION_TERMS];
} Derivatives;

int binomial(const int n, const int k);
float integer_powf(const float x, const unsigned int n);
void addPowers(thread Multipole &mp);

Multipole P2M(device int* treeStructure, device float* masses, device float3* positions, device float* grav, device bool* active, device int* nextActiveTime, device int* globalTime, int treePointer);
float3 L2P(Local local, float3 pos);
Multipole M2M(device int* treeStructure, device Multipole* multipoles, device bool* active, device uint* parentIndexes, uint index, int treePointer);
void L2L(device int* treeStructure, device Multipole* multipoles, device Local* locals, thread Multipole& mp, thread Local& local, int treePointer);
Derivatives derivatives(float3 vec, float eps);
Local M2L(float3 x, Multipole mp);
bool gravity_M2L_accept(Multipole A, Multipole B);




#endif /* poles_h */
