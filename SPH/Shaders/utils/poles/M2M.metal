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
