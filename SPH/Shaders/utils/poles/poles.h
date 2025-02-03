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

typedef struct {
    float mass;
    float3 pos;
    float3 min;
    float3 max;
    float size;
    float3 external;
} Multipole;

void addParticletoMonopole(thread Multipole &mp, float3 position, float mass);
void combineMultipole(thread Multipole &A, thread const Multipole &B);
void subtractMultipole(thread Multipole &A, thread const Multipole &B);
float3 multipoleAcceleration(thread const Multipole &mp, float3 r, float G);
void saveMultipole(device float* treeData, int index, Multipole mp);
void transformMultipole(thread Multipole &mp, thread const float3 &newR0);



//// A simple struct storing up to octupole moments, plus a reference-frame origin R0.
//typedef struct {
//    // Reference-frame origin for this multipole expansion:
//    float3 R0;
//
//    // Monopole:
//    float M;
//
//    // Dipole (3):
//    float3 D; // (Dx, Dy, Dz)
//
//    // Quadrupole (6 independent components, symmetric Q_{ij}):
//    float Qxx, Qxy, Qxz, Qyy, Qyz, Qzz;
//
//    // Octupole (10 independent components, fully symmetric O_{ijk}):
//    // We store them in a convenient order:
//    float Oxxx, Oxxy, Oxxz, Oxyy, Oxyz, Oxzz,
//          Oyyy, Oyyz, Oyzz, Ozzz;
//} Multipole;
//
//void addParticletoMonopole(thread Multipole &mp, float3 position, float mass);
//void combineMultipole(thread Multipole &A, thread const Multipole &B);
//void subtractMultipole(thread Multipole &A, thread const Multipole &B);
//float3 multipoleAcceleration(thread const Multipole &mp, float3 r, float G);
//Multipole getMultipole();
//Multipole getMultipole(device float* treeData, int index);
//void saveMultipole(device float* treeData, int index, Multipole mp);
//void transformMultipole(thread Multipole &mp, thread const float3 &newR0);

#endif /* poles_h */
