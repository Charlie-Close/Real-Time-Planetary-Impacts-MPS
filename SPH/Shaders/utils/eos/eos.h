//
//  eos.h
//  SPH
//
//  Created by Charlie Close on 03/11/2024.
//

#ifndef eos_h
#define eos_h

#include <metal_stdlib>
using namespace metal;

float pressureEos(float density, float internalEnergy, int materialId);

float speedOfSoundEos(float density, float internalEnergy, int materialId);

#endif /* eos_h */
