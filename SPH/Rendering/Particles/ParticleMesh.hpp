//
//  ParticleMesh.hpp
//  SPH
//
//  Created by Charlie Close on 24/01/2025.
//

#ifndef ParticleMesh_hpp
#define ParticleMesh_hpp

#include <stdio.h>
#include <vector>
#include <simd/simd.h>

std::tuple<std::vector<simd_float3>, std::vector<simd_float3>, std::vector<uint16_t>> generateSphere(float size, int subdivisions);

#endif /* ParticleMesh_hpp */
