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
    
    newMp.expansion[X] = 0;//mp.expansion[X] + r.x * mp.expansion[M];
    newMp.expansion[Y] = 0;//mp.expansion[Y] + r.y * mp.expansion[M];
    newMp.expansion[Z] = 0;//mp.expansion[Z] + r.z * mp.expansion[M];
    
#if P > 1
    newMp.expansion[XX] = mp.expansion[XX] + /*r.x * mp.expansion[X] +*/ 0.5 * r.x * r.x * mp.expansion[M];
    newMp.expansion[XY] = mp.expansion[XY] + /*r.y * mp.expansion[X] + r.x * mp.expansion[Y] +*/ r.x * r.y * mp.expansion[M];
    newMp.expansion[XZ] = mp.expansion[XZ] + /*r.z * mp.expansion[X] + r.x * mp.expansion[Z] +*/ r.x * r.z * mp.expansion[M];
    
    newMp.expansion[YY] = mp.expansion[YY] + /*r.y * mp.expansion[Y] +*/ 0.5 * r.y * r.y * mp.expansion[M];
    newMp.expansion[YZ] = mp.expansion[YZ] + /*r.z * mp.expansion[Y] + r.y * mp.expansion[Z] +*/ r.y * r.z * mp.expansion[M];
    
    newMp.expansion[ZZ] = mp.expansion[ZZ] + /*r.z * mp.expansion[Z] +*/ 0.5 * r.z * r.z * mp.expansion[Z];
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
#if P > 3
    /* Shift 4th order terms (1st order mpole (all 0) commented out) */
    newMp.expansion[ZZZZ] = mp.expansion[ZZZZ] + r.z * mp.expansion[ZZZ] +
                 0.5 * r.z * r.z * mp.expansion[ZZ] /* + 0.1666666666666667 * r.z * r.z * r.z * mp.expansion[Z] */ +
                 0.041666666666666667 * r.z * r.z * r.z * r.z * mp.expansion[M];
    newMp.expansion[YZZZ] = mp.expansion[YZZZ] + r.z * mp.expansion[YZZ] +
                 0.5 * r.z * r.z * mp.expansion[YZ] /* + 0.1666666666666667 * r.z * r.z * r.z * mp.expansion[Y] */ +
                 r.y * mp.expansion[ZZZ] +
                 r.y * r.z * mp.expansion[ZZ] /* + 0.5 * r.y * r.z * r.z * mp.expansion[Z] */ +
                 0.1666666666666667 * r.y * r.z * r.z * r.z * mp.expansion[M];
    newMp.expansion[YYZZ] = mp.expansion[YYZZ] + r.z * mp.expansion[YYZ] + 0.5 * r.z * r.z * mp.expansion[YY] +
                 r.y * mp.expansion[YZZ] +
                 r.y * r.z * mp.expansion[YZ] /* + 0.5 * r.y * r.z * r.z * mp.expansion[Y] */ +
                 0.5 * r.y * r.y * mp.expansion[ZZ] /* + 0.5 * r.y * r.y * r.z * mp.expansion[Z] */ +
                 0.25 * r.y * r.y * r.z * r.z * mp.expansion[M];
    newMp.expansion[YYYZ] = mp.expansion[YYZ] + r.z * mp.expansion[YYY] + r.y * mp.expansion[YYZ] +
                 r.y * r.z * mp.expansion[YY] +
                 0.5 * r.y * r.y * mp.expansion[YZ] /* + 0.5 * r.y * r.y * r.z * mp.expansion[Y] */
                                        /* + 0.1666666666666667 * r.y * r.y * r.y * mp.expansion[Z] */
                 + 0.1666666666666667 * r.y * r.y * r.y * r.z * mp.expansion[M];
    newMp.expansion[YYYY] = mp.expansion[YYYY] + r.y * mp.expansion[YYY] +
                 0.5 * r.y * r.y * mp.expansion[YY] /* + 0.1666666666666667 * r.y * r.y * r.y * mp.expansion[Y] */ +
                 0.041666666666666667 * r.y * r.y * r.y * r.y * mp.expansion[M];
    newMp.expansion[XZZZ] = mp.expansion[XZZZ] + r.z * mp.expansion[XZZ] +
                 0.5 * r.z * r.z * mp.expansion[XZ] /* + 0.1666666666666667 * r.z * r.z * r.z * mp.expansion[X] */ +
                 r.z * mp.expansion[ZZZ] +
                 r.x * r.z * mp.expansion[ZZ] /* + 0.5 * r.x * r.z * r.z * mp.expansion[Z] */ +
                 0.1666666666666667 * r.x * r.z * r.z * r.z * mp.expansion[M];
    newMp.expansion[XYZZ] = mp.expansion[XYZZ] + r.z * mp.expansion[XYZ] + 0.5 * r.z * r.z * mp.expansion[XY] +
                 r.y * mp.expansion[XZZ] +
                 r.y * r.z * mp.expansion[XZ] /* + 0.5 * r.y * r.z * r.z * mp.expansion[X] */ +
                 r.z * mp.expansion[YZZ] +
                 r.x * r.z * mp.expansion[YZ] /* + 0.5 * r.x * r.z * r.z * mp.expansion[Y] */ +
                 r.x * r.y * mp.expansion[ZZ] /* + r.x * r.y * r.z * mp.expansion[Z] */ +
                 0.5 * r.x * r.y * r.z * r.z * mp.expansion[M];
    newMp.expansion[XYYZ] = mp.expansion[XYYZ] + r.z * mp.expansion[XYY] + r.y * mp.expansion[XYZ] +
                 r.y * r.z * mp.expansion[XY] +
                 0.5 * r.y * r.y * mp.expansion[XZ] /* + 0.5 * r.y * r.y * r.z * mp.expansion[X] */ +
                 r.z * mp.expansion[YYZ] + r.x * r.z * mp.expansion[YY] +
                 r.x * r.y * mp.expansion[YZ] /* + r.x * r.y * r.z * mp.expansion[Y] */
                                        /* + 0.5 * r.x * r.y * r.y * mp.expansion[Z] */
                 + 0.5 * r.x * r.y * r.y * r.z * mp.expansion[M];
    newMp.expansion[XYYY] = mp.expansion[XYYY] + r.y * mp.expansion[XYY] +
                 0.5 * r.y * r.y * mp.expansion[XY] /* + 0.1666666666666667 * r.y * r.y * r.y * mp.expansion[X] */ +
                 r.z * mp.expansion[YYY] +
                 r.x * r.y * mp.expansion[YY] /* + 0.5 * r.x * r.y * r.y * mp.expansion[Y] */ +
                 0.1666666666666667 * r.x * r.y * r.y * r.y * mp.expansion[M];
    newMp.expansion[XXZZ] = mp.expansion[XXZZ] + r.z * mp.expansion[XXZ] + 0.5 * r.z * r.z * mp.expansion[XX] +
                 r.z * mp.expansion[XZZ] +
                 r.x * r.z * mp.expansion[XZ] /* + 0.5 * r.x * r.z * r.z * mp.expansion[X] */ +
                 0.5 * r.x * r.x * mp.expansion[ZZ] /* + 0.5 * r.x * r.x * r.z * mp.expansion[Z] */ +
                 0.25 * r.x * r.x * r.z * r.z * mp.expansion[M];
    newMp.expansion[XXYZ] = mp.expansion[XXYZ] + r.z * mp.expansion[XXY] + r.y * mp.expansion[XXZ] +
                 r.y * r.z * mp.expansion[XX] + r.z * mp.expansion[XYZ] +
                 r.x * r.z * mp.expansion[XY] +
                 r.x * r.y * mp.expansion[XZ] /* + r.x * r.y * r.z * mp.expansion[X] */ +
                 0.5 * r.x * r.x * mp.expansion[YZ] /* + 0.5 * r.x * r.x * r.z * mp.expansion[Y] */
                                        /* + 0.5 * r.x * r.x * r.y * mp.expansion[Z] */
                 + 0.5 * r.x * r.x * r.y * r.z * mp.expansion[M];
    newMp.expansion[XXYY] = mp.expansion[XXYY] + r.y * mp.expansion[XXY] + 0.5 * r.y * r.y * mp.expansion[XX] +
                 r.z * mp.expansion[XYY] +
                 r.x * r.y * mp.expansion[XY] /* + 0.5 * r.x * r.y * r.y * mp.expansion[X] */ +
                 0.5 * r.x * r.x * mp.expansion[YY] /* + 0.5 * r.x * r.x * r.y * mp.expansion[Y] */ +
                 0.25 * r.x * r.x * r.y * r.y * mp.expansion[M];
    newMp.expansion[XXXZ] = mp.expansion[XXXZ] + r.z * mp.expansion[XXX] + r.z * mp.expansion[XXZ] +
                 r.x * r.z * mp.expansion[XX] +
                 0.5 * r.x * r.x * mp.expansion[XZ] /* + 0.5 * r.x * r.x * r.z * mp.expansion[X] */
                                        /* + 0.1666666666666667 * r.x * r.x * r.x * mp.expansion[Z] */
                 + 0.1666666666666667 * r.x * r.x * r.x * r.z * mp.expansion[M];
    newMp.expansion[XXXY] = mp.expansion[XXXY] + r.y * mp.expansion[XXX] + r.z * mp.expansion[XXY] +
                 r.x * r.y * mp.expansion[XX] +
                 0.5 * r.x * r.x * mp.expansion[XY] /* + 0.5 * r.x * r.x * r.y * mp.expansion[X] */
                                        /* + 0.1666666666666667 * r.x * r.x * r.x * mp.expansion[Y] */
                 + 0.1666666666666667 * r.x * r.x * r.x * r.y * mp.expansion[M];
    newMp.expansion[XXXX] = mp.expansion[XXXX] + r.z * mp.expansion[XXX] +
                 0.5 * r.x * r.x * mp.expansion[XX] /* + 0.1666666666666667 * r.x * r.x * r.x * mp.expansion[X] */ +
                 0.041666666666666667 * r.x * r.x * r.x * r.x * mp.expansion[M];
#endif
    return newMp;
}

Multipole M2M(device int* treeStructure, device Multipole* multipoles, device bool* active, device unsigned long* parentIndexes, uint index, int treePointer) {
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
    mp.max = float3(-MAXFLOAT);
    mp.min = float3(MAXFLOAT);
    mp.eta = 0;
    
    // Get the start and end of our child pointers (8 of them starting 2 ahead of treePointer)
    int start = treePointer + 2;
    int end = start + 8;
    bool allLeaves = true;
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
        mp.eta = max(mp.eta, childMp.eta);
        // If at least one child is active, we are active (used to skip inactive
        // nodes in down pass)
        if (!mp.active and childMp.active) {
            mp.active = true;
        }
        if (allLeaves and treeStructure[childPointer] == 0) {
            allLeaves = false;
        }
        
        mp.max = max(mp.max, childMp.max);
        mp.min = min(mp.min, childMp.min);
    }
    
    float3 dims = mp.max - mp.min;
    mp.size = max3(dims.x, dims.y, dims.z);
    mp.pos /= mass;
    
    // Second pass: get our multipole expansion. If all are children are leaves, we sync them
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
        
        if (allLeaves and mp.active and not childMp.active) {
            // Activate child to sync it with neighbouring leaves. Improves time stepping
            multipoles[childDataPointer].active = true;
            int nParticles = treeStructure[childPointer];
            int childStart = childPointer + 2;
            int childEnd = start + nParticles;
            for (int j = childStart; j < childEnd; j++) {
                active[treeStructure[j]] = true;
            }
        }
    }
    // This is for multipole acceptance criterion
    addPowers(mp);
    
    return mp;
}
