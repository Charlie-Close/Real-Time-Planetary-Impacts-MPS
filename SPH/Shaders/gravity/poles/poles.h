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

constant uint XXXX = 20;
constant uint XXXY = 21;
constant uint XXXZ = 22;
constant uint XXYY = 23;
constant uint XXYZ = 24;
constant uint XXZZ = 25;
constant uint XYYY = 26;
constant uint XYYZ = 27;
constant uint XYZZ = 28;
constant uint XZZZ = 29;
constant uint YYYY = 30;
constant uint YYYZ = 31;
constant uint YYZZ = 32;
constant uint YZZZ = 33;
constant uint ZZZZ = 34;

typedef struct {
    float3 pos;
    float3 min;
    float3 max;
    float size;
    float expansion[N_EXPANSION_TERMS];
    float power[P+1];
    float minGrav;
    bool active;
    float eta;
} Multipole;


typedef struct {
    float3 pos;
    float expansion[N_EXPANSION_TERMS];
} Local;

typedef struct {
    float expansion[N_EXPANSION_TERMS];
} Derivatives;

float integer_powf(const float x, const unsigned int n);
void addPowers(thread Multipole &mp);

Multipole P2M(device int* treeStructure, device float* masses, device float3* positions, device float* grav, device float* h, device bool* active, device int* nextActiveTime, device int* globalTime, device int& dt, int treePointer);
float3 L2P(Local local, float3 pos);
Multipole M2M(device int* treeStructure, device Multipole* multipoles, device bool* active, device unsigned long* parentIndexes, uint index, int treePointer);
void L2L(device int* treeStructure, device Multipole* multipoles, device Local* locals, thread Multipole& mp, thread Local& local, int treePointer);
Derivatives derivatives(float3 vec, float eps);
Local M2L(float3 x, Multipole mp);
bool gravity_M2L_accept(Multipole A, Multipole B);

constant int binomial_coeffs[9][9] = {
    {1, 0, 0, 0, 0, 0, 0, 0, 0},     {1, 1, 0, 0, 0, 0, 0, 0, 0},
    {1, 2, 1, 0, 0, 0, 0, 0, 0},     {1, 3, 3, 1, 0, 0, 0, 0, 0},
    {1, 4, 6, 4, 1, 0, 0, 0, 0},     {1, 5, 10, 10, 5, 1, 0, 0, 0},
    {1, 6, 15, 20, 15, 6, 1, 0, 0},  {1, 7, 21, 35, 35, 21, 7, 1, 0},
    {1, 8, 28, 56, 70, 56, 28, 8, 1}

};




#endif /* poles_h */
