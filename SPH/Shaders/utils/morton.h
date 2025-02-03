//
//  morton.h
//  SPH
//
//  Created by Charlie Close on 02/02/2025.
//

#ifndef morton_h
#define morton_h

#include <metal_stdlib>
using namespace metal;

// -------------- GPU-friendly bit ops ---------------
// expandBits_32: expand a 10-bit input into 30 bits with 2 zeroes between each bit.
static inline uint expandBits_32(uint v)
{
    v = (v * 0x00010001u) & 0xFF0000FFu; // spread lower 8 bits into pairs of zeros
    v = (v * 0x00000101u) & 0x0F00F00Fu; // spread
    v = (v * 0x00000011u) & 0xC30C30C3u; // spread
    v = (v * 0x00000005u) & 0x49249249u; // spread
    return v;
}

// morton3D: interleave bits of x, y, z
static inline uint morton3D(uint x, uint y, uint z)
{
    // Expand each into 30 bits with zero spacing
    uint xx = expandBits_32(x);
    uint yy = expandBits_32(y) << 1; // shift y bits one to the left
    uint zz = expandBits_32(z) << 2; // shift z bits two to the left
    // Combine
    return (xx | yy | zz);
}

#endif /* morton_h */
