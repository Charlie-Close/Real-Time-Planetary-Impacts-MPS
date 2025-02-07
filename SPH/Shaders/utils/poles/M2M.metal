//
//  M2M.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include "poles.h"
#include <metal_stdlib>
using namespace metal;

Multipole transformMultipole(Multipole mp, float3 r) {
    Multipole newMp;
    newMp.pos = mp.pos + r;
    newMp.expansion[M] = mp.expansion[M];
    
    newMp.expansion[X] = mp.expansion[X] + r.x * mp.expansion[M];
    newMp.expansion[Y] = mp.expansion[Y] + r.y * mp.expansion[M];
    newMp.expansion[Z] = mp.expansion[Z] + r.z * mp.expansion[M];
    
#if P > 1
    newMp.expansion[XX] = mp.expansion[XX] + r.x * mp.expansion[X] + 0.5 * r.x * r.x * mp.expansion[M];
    newMp.expansion[XY] = mp.expansion[XY] + r.y * mp.expansion[X] + r.x * mp.expansion[Y] + r.x * r.y * mp.expansion[M];
    newMp.expansion[XZ] = mp.expansion[XZ] + r.z * mp.expansion[X] + r.x * mp.expansion[Z] + r.x * r.z * mp.expansion[M];
    
    newMp.expansion[YY] = mp.expansion[YY] + r.y * mp.expansion[Y] + 0.5 * r.y * r.y * mp.expansion[M];
    newMp.expansion[YZ] = mp.expansion[YZ] + r.z * mp.expansion[Y] + r.y * mp.expansion[Z] + r.y * r.z * mp.expansion[M];
    
    newMp.expansion[ZZ] = mp.expansion[ZZ] + r.z * mp.expansion[Z] + 0.5 * r.z * r.z * mp.expansion[Z];
#endif
#if P > 2
    /* Shift 3rd order terms (1st order mpole (all 0) commented out) */
    newMp.expansion[ZZZ] = mp.expansion[ZZZ] +
                 r.z * mp.expansion[ZZ] /* + 0.5 * r.z * r.z * mp.expansion[];//M_001 */ +
                 0.166666666667 * r.z * r.z * r.z * mp.expansion[M];
    newMp.expansion[YZZ] = mp.expansion[YZZ] +
                 r.z * mp.expansion[YZ] /* + 0.5 * r.z * r.z * mp.expansion[];//M_010 */ +
                 r.y * mp.expansion[ZZ] /* + r.y * r.z * mp.expansion[];//M_001 */ +
                 0.5 * r.y * r.z * r.z * mp.expansion[M];
    newMp.expansion[YYZ] = mp.expansion[YYZ] + r.z * mp.expansion[YY] +
                 r.y * mp.expansion[YZ] /* + r.y * r.z * mp.expansion[];//M_010 */
                                        /* + 0.5 * r.y * r.y * mp.expansion[];//M_001 */
                 + 0.5 * r.y * r.y * r.z * mp.expansion[M];
    newMp.expansion[YYY] = mp.expansion[YYY] +
                 r.y * mp.expansion[YY] /* + 0.5 * r.y * r.y * mp.expansion[];//M_010 */ +
                 0.166666666667 * r.y * r.y * r.y * mp.expansion[M];
    newMp.expansion[XZZ] = mp.expansion[XZZ] +
                 r.z * mp.expansion[XZ] /* + 0.5 * r.z * r.z * mp.expansion[];//M_100 */ +
                 r.x * mp.expansion[ZZ] /* + r.x * r.z * mp.expansion[];//M_001 */ +
                 0.5 * r.x * r.z * r.z * mp.expansion[M];
    newMp.expansion[XYZ] = mp.expansion[XYZ] + r.z * mp.expansion[XY] +
                 r.y * mp.expansion[XZ] /* + r.y * r.z * mp.expansion[];//M_100 */ +
                 r.x * mp.expansion[YZ] /* + r.x * r.z * mp.expansion[];//M_010 */
                                        /* + r.x * r.y * mp.expansion[];//M_001 */
                 + r.x * r.y * r.z * mp.expansion[M];
    newMp.expansion[XYY] = mp.expansion[XYY] +
                 r.y * mp.expansion[XY] /* + 0.5 * r.y * r.y * mp.expansion[];//M_100 */ +
                 r.x * mp.expansion[YY] /* + r.x * r.y * mp.expansion[];//M_010 */ +
                 0.5 * r.x * r.y * r.y * mp.expansion[M];
    newMp.expansion[XXZ] = mp.expansion[XXZ] + r.z * mp.expansion[XX] +
                 r.x * mp.expansion[XZ] /* + r.x * r.z * mp.expansion[];//M_100 */
                                        /* + 0.5 * r.x * r.x * mp.expansion[];//M_001 */
                 + 0.5 * r.x * r.x * r.z * mp.expansion[M];
    newMp.expansion[XXY] = mp.expansion[XXY] + r.y * mp.expansion[XX] +
                 r.x * mp.expansion[XY] /* + r.x * r.y * mp.expansion[];//M_100 */
                                        /* + 0.5 * r.x * r.x * mp.expansion[];//M_010 */
                 + 0.5 * r.x * r.x * r.y * mp.expansion[M];
    newMp.expansion[XXX] = mp.expansion[XXX] +
                 r.x * mp.expansion[XX] /* + 0.5 * r.x * r.x * mp.expansion[];//M_100 */ +
                 0.1666666667 * r.x * r.x * r.x * mp.expansion[M];
#endif
    return newMp;
}

Multipole M2M(device int* treeStructure, device Multipole* multipoles, device uint* parentIndexes, uint index, int treePointer) {
    // If we are looking at a branch (signaled by nParticles == 0), we go through all
    // 8 child nodes twice. The first time to find our center of mass, and the second
    // to sum our children's multipoles
    
    // Initialise our multipole with everything 0, and minGrav as MAXFLOAT
    Multipole mp;
    mp.minGrav = MAXFLOAT;
    mp.pos = { 0, 0, 0 };
    mp.active = false;
    for (uint i = 0; i < N_EXPANSION_TERMS; i++) {
        mp.expansion[i] = 0;
    }
    float mass = 0;
    bool first = true;
    
    // Get the start and end of our child pointers (8 of them starting 2 ahead of treePointer)
    int start = treePointer + 2;
    int end = start + 8;
    // First pass: we find COM, min and max coords and whether we are active
    for (int i = start; i < end; i++) {
        int childPointer = treeStructure[i];
        // -1 signifies there is not a child here so we skip
        if (childPointer == -1) {
            continue;
        }
        int childDataPointer = treeStructure[childPointer + 1];
        // We let the child know where to find us (used in down pass)
        parentIndexes[childDataPointer] = index;
        Multipole childMp = multipoles[childDataPointer];
        // Add mass and position * mass so we can calculate our COM
        mp.minGrav = min(childMp.minGrav, mp.minGrav);
        mass += childMp.expansion[M];
        mp.pos += childMp.pos * childMp.expansion[M];
        // If at least one child is active, we are active (used to skip inactive
        // nodes in down pass)
        if (!mp.active && childMp.active) {
            mp.active = true;
        }
        
        // This is kinda jank, but you get the idea - find min and max coords
        if (first) {
            mp.max = childMp.max;
            mp.min = childMp.min;
            first = false;
        } else {
            mp.max.x = max(childMp.max.x, mp.max.x);
            mp.max.y = max(childMp.max.y, mp.max.y);
            mp.max.z = max(childMp.max.z, mp.max.z);
            mp.min.x = min(childMp.min.x, mp.min.x);
            mp.min.y = min(childMp.min.y, mp.min.y);
            mp.min.z = min(childMp.min.z, mp.min.z);
        }
    }
    
    float3 dims = mp.max - mp.min;
    mp.size = max3(dims.x, dims.y, dims.z);
    mp.pos /= mass;
    
    // Second pass: get our multipole expansion
    for (int i = start; i < end; i++) {
        int childPointer = treeStructure[i];
        // -1 signifies there is not a child here so we skip
        if (childPointer == -1) {
            continue;
        }
        int childDataPointer = treeStructure[childPointer + 1];
        Multipole childMp = multipoles[childDataPointer];
        // Now we know our COM, we can shift the child's multipole expansion
        // to center at our COM, and then we can just sum the expansion terms
        float3 r = mp.pos - childMp.pos;
        Multipole transformed = transformMultipole(childMp, r);
        for (uint j = 0; j < N_EXPANSION_TERMS; j++) {
            mp.expansion[j] += transformed.expansion[j];
        }
    }
    // This is for multipole acceptance criterion
    addPowers(mp);
    
    return mp;
}
