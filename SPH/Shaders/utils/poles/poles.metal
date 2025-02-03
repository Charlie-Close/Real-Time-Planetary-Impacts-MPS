//
//  poles.metal
//  SPH
//
//  Created by Charlie Close on 23/01/2025.
//

#include "poles.h"





//Multipole getMultipole() {
//    Multipole mp;
//    mp.R0 = { 0, 0, 0 };
//    mp.M = 0;
//    mp.D = { 0, 0, 0 };
//    mp.Qxx = mp.Qxy = mp.Qxz = mp.Qyy = mp.Qyz = mp.Qzz = 0.f;
//    mp.Oxxx = mp.Oxxy = mp.Oxxz = mp.Oxyy = mp.Oxyz = mp.Oxzz = 0.f;
//    mp.Oyyy = mp.Oyyz = mp.Oyzz = mp.Ozzz = 0.f;
//    
//    return mp;
//};
//
//Multipole getMultipole(device float* treeData, int index) {
//    Multipole mp;
//    mp.M = treeData[index];
//    mp.D = { treeData[index + 1], treeData[index + 2], treeData[index + 3] };
//    mp.Qxx = treeData[index + 4]; mp.Qxy = treeData[index + 5]; mp.Qxz = treeData[index + 6]; mp.Qyy = treeData[index + 7]; mp.Qyz = treeData[index + 8]; mp.Qzz = treeData[index + 9];
//    mp.Oxxx = treeData[index + 10]; mp.Oxxy = treeData[index + 11]; mp.Oxxz = treeData[index + 12]; mp.Oxyy = treeData[index + 13]; mp.Oxyz = treeData[index + 14]; mp.Oxzz = treeData[index + 15]; mp.Oyyy = treeData[index + 16]; mp.Oyyz = treeData[index + 17]; mp.Oyzz = treeData[index + 18]; mp.Ozzz = treeData[index + 19];
//    mp.R0 = { treeData[index + 20], treeData[index + 21], treeData[index + 22] };
//
//    return mp;
//};
//
//void saveMultipole(device float* treeData, int index, Multipole mp) {
//    treeData[index] = mp.M;
//    treeData[index + 1] = mp.D[0]; treeData[index + 2] = mp.D[1]; treeData[index + 3] = mp.D[2];
//    treeData[index + 4] = mp.Qxx; treeData[index + 5] = mp.Qxy; treeData[index + 6] = mp.Qxz; treeData[index + 7] = mp.Qyy; treeData[index + 8] = mp.Qyz; treeData[index + 9] = mp.Qzz;
//    treeData[index + 10] = mp.Oxxx; treeData[index + 11] = mp.Oxxy; treeData[index + 12] = mp.Oxxz; treeData[index + 13] = mp.Oxyy; treeData[index + 14] = mp.Oxyz; treeData[index + 15] = mp.Oxzz; treeData[index + 16] = mp.Oyyy; treeData[index + 17] = mp.Oyyz; treeData[index + 18] = mp.Oyzz; treeData[index + 19] = mp.Ozzz;
//    treeData[index + 20] = mp.R0[0]; treeData[index + 21] = mp.R0[1]; treeData[index + 22] = mp.R0[2];
//
//};
//
//void addParticletoMonopole(thread Multipole &mp, float3 position, float mass) {
//    // Relative position w.r.t. the center of mass
//    float3 dr = position - mp.R0;
//    
//    // Dipole
//    mp.D += mass * dr;
//
//    // Quadrupole (symmetric matrix)
//    mp.Qxx += mass * dr.x * dr.x;
//    mp.Qxy += mass * dr.x * dr.y;
//    mp.Qxz += mass * dr.x * dr.z;
//    mp.Qyy += mass * dr.y * dr.y;
//    mp.Qyz += mass * dr.y * dr.z;
//    mp.Qzz += mass * dr.z * dr.z;
//
//    // Octupole (fully symmetric tensor)
//    mp.Oxxx += mass * dr.x * dr.x * dr.x;
//    mp.Oxxy += mass * dr.x * dr.x * dr.y;
//    mp.Oxxz += mass * dr.x * dr.x * dr.z;
//    mp.Oxyy += mass * dr.x * dr.y * dr.y;
//    mp.Oxyz += mass * dr.x * dr.y * dr.z;
//    mp.Oxzz += mass * dr.x * dr.z * dr.z;
//    mp.Oyyy += mass * dr.y * dr.y * dr.y;
//    mp.Oyyz += mass * dr.y * dr.y * dr.z;
//    mp.Oyzz += mass * dr.y * dr.z * dr.z;
//    mp.Ozzz += mass * dr.z * dr.z * dr.z;
//}
//
//// Helper to combine B into A, using A's reference frame.
//// That is: shift B's multipole expansion from R0^B to R0^A, then add to A.
//void combineMultipole(thread Multipole &A, thread const Multipole &B)
//{
//    // Compute shift vector Δ = (B.R0 - A.R0).
//    float3 dR = (B.R0 - A.R0);
//
//    //---------------------------------------------------------
//    // 1) Monopole
//    //---------------------------------------------------------
//    float Mtemp = A.M + B.M;
//
//    //---------------------------------------------------------
//    // 2) Dipole
//    //---------------------------------------------------------
//    // SHIFTED dipole from B:
//    //    D_B_shift = D_B + M_B * dR
//    float3 D_B_shifted = B.D + (B.M * dR);
//
//    // Add to A:
//    float3 Dtemp = A.D + D_B_shifted;
//
//    //---------------------------------------------------------
//    // 3) Quadrupole
//    //---------------------------------------------------------
//    // We'll define a small helper for each component:
//    // Q'_xx = QA_xx + QB_xx + ( Δ.x*DB.x + Δ.x*DB.x ) + M_B*Δ.x^2 ...
//    // But more systematically:
//    //   QA_xx + [QB_xx + Δ.x*DB.x + Δ.x*DB.x + M_B*Δ.x^2]
//    //   => +2 * Δ.x * D_B.x ...
//    // Similarly for xy, xz, yy, yz, zz.
//
//    float Qxx_temp = A.Qxx + B.Qxx
//                     + 2.f * dR.x * B.D.x + B.M * (dR.x * dR.x);
//    float Qxy_temp = A.Qxy + B.Qxy
//                     + dR.x * B.D.y + dR.y * B.D.x
//                     + B.M * (dR.x * dR.y);
//    float Qxz_temp = A.Qxz + B.Qxz
//                     + dR.x * B.D.z + dR.z * B.D.x
//                     + B.M * (dR.x * dR.z);
//    float Qyy_temp = A.Qyy + B.Qyy
//                     + 2.f * dR.y * B.D.y + B.M * (dR.y * dR.y);
//    float Qyz_temp = A.Qyz + B.Qyz
//                     + dR.y * B.D.z + dR.z * B.D.y
//                     + B.M * (dR.y * dR.z);
//    float Qzz_temp = A.Qzz + B.Qzz
//                     + 2.f * dR.z * B.D.z + B.M * (dR.z * dR.z);
//
//    //---------------------------------------------------------
//    // 4) Octupole
//    //---------------------------------------------------------
//    // We apply:
//    //  O'_{ijk} = O^A_{ijk} +
//    //             [ O^B_{ijk}
//    //               + Δ_i Q^B_{jk} + Δ_j Q^B_{ik} + Δ_k Q^B_{ij}
//    //               + Δ_i Δ_j D^B_k + Δ_j Δ_k D^B_i + Δ_k Δ_i D^B_j
//    //               + M_B Δ_i Δ_j Δ_k ]
//    //
//    // We'll do it component by component.
//    //
//    // For example, Oxxx' = A.Oxxx + [ B.Oxxx
//    //    + dR.x * B.Qxx + dR.x * B.Qxx + dR.x * B.Qxx  (the 3 ways i=j=k=x)
//    //    + dR.x * dR.x * B.D.x + ...
//    //    + B.M * dR.x^3
//    // ]
//    // Actually we have to be careful to pick out the right Q_{jk}, etc.
//
//    auto Q_B = [&](float3 r)->float3 {
//        // For "Q*r":
//        //   (Qxx*r.x + Qxy*r.y + Qxz*r.z,
//        //    Qxy*r.x + Qyy*r.y + Qyz*r.z,
//        //    Qxz*r.x + Qyz*r.y + Qzz*r.z)
//        return float3{ B.Qxx*r.x + B.Qxy*r.y + B.Qxz*r.z,
//                       B.Qxy*r.x + B.Qyy*r.y + B.Qyz*r.z,
//                       B.Qxz*r.x + B.Qyz*r.y + B.Qzz*r.z };
//    };
//
//    // For each octupole component, we add:
//    //   B.Oxxx, Oxxy, ...
//    // plus the terms with Q_B + D_B, etc.
//    // We'll define a small function "Q_ijk_shift" that calculates
//    // (Δ_i Q^B_jk + Δ_j Q^B_ik + Δ_k Q^B_ij).
//
//    // Let's do it explicitly to avoid confusion:
//
//    // 1) Oxxx
//    float Oxxx_temp = A.Oxxx + B.Oxxx;
//
//    // Add Q shifts => we want i=j=k=x:
//    //   Δ.x * Q_{xx}^B + Δ.x * Q_{xx}^B + Δ.x * Q_{xx}^B = 3 * Δ.x * Q_{xx}^B
//    // But we have Qxx, Qxy,... so Q_{xx}^B = B.Qxx
//    Oxxx_temp += 3.f * (dR.x * B.Qxx);
//
//    // Add D shifts => i=j => x, x => Δ.x * Δ.x * D_B.x => Δ.x^2 * B.D.x
//    // times 3 because (i=j=k=x => permutations), effectively 1 combination for (x,x,x).
//    Oxxx_temp += 3.f * (dR.x*dR.x * B.D.x);
//
//    // Add M shifts => M_B * Δ.x^3
//    Oxxx_temp += B.M * (dR.x * dR.x * dR.x);
//
//    // 2) Oxxy
//    float Oxxy_temp = A.Oxxy + B.Oxxy;
//
//    // Q shift => i=j=x, k=y =>
//    //   + Δ.x * Q_xy^B   (i=x => Δ.x * Q_{xy})
//    //   + Δ.x * Q_xy^B   (j=x => Δ.x * Q_{xy})
//    //   + Δ.y * Q_{xx}^B (k=y => Δ.y * Q_{xx})
//    // So total: Q_{xy}^B*(Δ.x+Δ.x) + Q_{xx}^B*Δ.y => 2 Δ.x Q_{xy}^B + Δ.y Q_{xx}^B
//    Oxxy_temp += 2.f*dR.x*B.Qxy + dR.y*B.Qxx;
//
//    // D shift =>
//    //   (Δ.x Δ.x D_B.y) + (Δ.x Δ.y D_B.x) + (Δ.y Δ.x D_B.x)
//    // Because i=j=x => that covers Δ.x Δ.x D_B.y, j=x,k=y => etc.
//    // Carefully enumerating:
//    //   i=x, j=x => we add Δ.x*Δ.x * D_B.(k)?  k=y => D_B.y
//    //   i=x, k=y => Δ.x*Δ.y * D_B.x
//    //   j=x, k=y => Δ.x*Δ.y * D_B.x
//    // Summation => Δ.x^2 * D_B.y + 2*(Δ.x*Δ.y*D_B.x)
//    Oxxy_temp += (dR.x*dR.x * B.D.y) + 2.f*(dR.x*dR.y * B.D.x);
//
//    // M shift => M_B Δ.x^2 Δ.y
//    Oxxy_temp += B.M * (dR.x*dR.x*dR.y);
//
//    // ... Similarly for the other 8 octupole components.
//    // Rather than type them all out in text, let's just show them in code below:
//
//    float Oxxz_temp = A.Oxxz + B.Oxxz;
//    {
//        // i=j=x, k=z
//        // Q shift => 2*Δ.x*Q_xz + Δ.z*Q_xx
//        Oxxz_temp += 2.f*dR.x*B.Qxz + dR.z*B.Qxx;
//        // D shift =>
//        //   i=x,j=x => Δ.x^2 * D_B.z
//        //   i=x,k=z => Δ.x*Δ.z * D_B.x
//        //   j=x,k=z => Δ.x*Δ.z * D_B.x
//        Oxxz_temp += (dR.x*dR.x * B.D.z) + 2.f*(dR.x*dR.z * B.D.x);
//        // M shift
//        Oxxz_temp += B.M*(dR.x*dR.x*dR.z);
//    }
//
//    float Oxyy_temp = A.Oxyy + B.Oxyy;
//    {
//        // i=x, j=y, k=y
//        // Q shift =>
//        //   i=x => Δ.x * Q_{yy}
//        //   j=y => Δ.y * Q_{xy}
//        //   k=y => Δ.y * Q_{xy}
//        // total => Δ.x*B.Qyy + 2.f*(Δ.y*B.Qxy)
//        Oxyy_temp += dR.x*B.Qyy + 2.f*(dR.y*B.Qxy);
//
//        // D shift =>
//        //   i=x,j=y => Δ.x*Δ.y * D_B.y
//        //   i=x,k=y => Δ.x*Δ.y * D_B.y
//        //   j=y,k=y => Δ.y^2 * D_B.x
//        // total => 2*(dR.x*dR.y*B.D.y) + dR.y*dR.y*B.D.x
//        Oxyy_temp += 2.f*(dR.x*dR.y*B.D.y) + (dR.y*dR.y*B.D.x);
//
//        // M shift => M_B * (dR.x*dR.y^2)
//        Oxyy_temp += B.M*(dR.x*dR.y*dR.y);
//    }
//
//    float Oxyz_temp = A.Oxyz + B.Oxyz;
//    {
//        // i=x, j=y, k=z
//        // Q shift =>
//        //   i=x => Δ.x * Q_{yz}
//        //   j=y => Δ.y * Q_{xz}
//        //   k=z => Δ.z * Q_{xy}
//        Oxyz_temp += dR.x*B.Qyz + dR.y*B.Qxz + dR.z*B.Qxy;
//        // D shift =>
//        //   i=x,j=y => Δ.x*Δ.y * D_B.z
//        //   i=x,k=z => Δ.x*Δ.z * D_B.y
//        //   j=y,k=z => Δ.y*Δ.z * D_B.x
//        Oxyz_temp += (dR.x*dR.y*B.D.z)
//                   + (dR.x*dR.z*B.D.y)
//                   + (dR.y*dR.z*B.D.x);
//        // M shift => M_B dR.x*dR.y*dR.z
//        Oxyz_temp += B.M*(dR.x*dR.y*dR.z);
//    }
//
//    float Oxzz_temp = A.Oxzz + B.Oxzz;
//    {
//        // i=x, j=z, k=z
//        // Q shift => Δ.x*Q_zz + 2.f*(Δ.z*Q_xz)
//        Oxzz_temp += dR.x*B.Qzz + 2.f*(dR.z*B.Qxz);
//        // D shift =>
//        //   i=x,j=z => Δ.x*Δ.z * D_B.z
//        //   i=x,k=z => Δ.x*Δ.z * D_B.z
//        //   j=z,k=z => Δ.z^2 * D_B.x
//        Oxzz_temp += 2.f*(dR.x*dR.z*B.D.z) + (dR.z*dR.z*B.D.x);
//        // M shift => B.M * (dR.x*dR.z^2)
//        Oxzz_temp += B.M*(dR.x*dR.z*dR.z);
//    }
//
//    float Oyyy_temp = A.Oyyy + B.Oyyy;
//    {
//        // i=j=k=y
//        // Q shift => 3.f * Δ.y * Q_{yy}
//        Oyyy_temp += 3.f*(dR.y*B.Qyy);
//        // D shift => 3.f*(dR.y*dR.y * B.D.y)
//        Oyyy_temp += 3.f*(dR.y*dR.y * B.D.y);
//        // M shift => B.M * dR.y^3
//        Oyyy_temp += B.M*(dR.y*dR.y*dR.y);
//    }
//
//    float Oyyz_temp = A.Oyyz + B.Oyyz;
//    {
//        // i=j=y, k=z
//        // Q shift => 2.f*(dR.y*B.Qyz) + dR.z*B.Qyy
//        Oyyz_temp += 2.f*(dR.y*B.Qyz) + dR.z*B.Qyy;
//        // D shift =>
//        //   i=y,j=y => dR.y^2 * D_B.z
//        //   i=y,k=z => dR.y*dR.z * D_B.y
//        //   j=y,k=z => dR.y*dR.z * D_B.y
//        Oyyz_temp += (dR.y*dR.y*B.D.z) + 2.f*(dR.y*dR.z*B.D.y);
//        // M shift => B.M * (dR.y^2*dR.z)
//        Oyyz_temp += B.M*(dR.y*dR.y*dR.z);
//    }
//
//    float Oyzz_temp = A.Oyzz + B.Oyzz;
//    {
//        // i=y, j=z, k=z
//        // Q shift => dR.y*B.Qzz + 2.f*(dR.z*B.Qyz)
//        Oyzz_temp += dR.y*B.Qzz + 2.f*(dR.z*B.Qyz);
//        // D shift =>
//        //   i=y,j=z => dR.y*dR.z * D_B.z
//        //   i=y,k=z => dR.y*dR.z * D_B.z
//        //   j=z,k=z => dR.z^2 * D_B.y
//        Oyzz_temp += 2.f*(dR.y*dR.z*B.D.z) + (dR.z*dR.z*B.D.y);
//        // M shift => B.M * (dR.y*dR.z^2)
//        Oyzz_temp += B.M*(dR.y*dR.z*dR.z);
//    }
//
//    float Ozzz_temp = A.Ozzz + B.Ozzz;
//    {
//        // i=j=k=z
//        // Q shift => 3.f * dR.z * B.Qzz
//        Ozzz_temp += 3.f*(dR.z*B.Qzz);
//        // D shift => 3.f*(dR.z*dR.z * B.D.z)
//        Ozzz_temp += 3.f*(dR.z*dR.z * B.D.z);
//        // M shift => B.M*(dR.z^3)
//        Ozzz_temp += B.M*(dR.z*dR.z*dR.z);
//    }
//
//    //---------------------------------------------------------
//    // Write them back into A.
//    //---------------------------------------------------------
//    A.M = Mtemp;
//    A.D = Dtemp;
//
//    A.Qxx = Qxx_temp;
//    A.Qxy = Qxy_temp;
//    A.Qxz = Qxz_temp;
//    A.Qyy = Qyy_temp;
//    A.Qyz = Qyz_temp;
//    A.Qzz = Qzz_temp;
//
//    A.Oxxx = Oxxx_temp;  A.Oxxy = Oxxy_temp;  A.Oxxz = Oxxz_temp;
//    A.Oxyy = Oxyy_temp;  A.Oxyz = Oxyz_temp;  A.Oxzz = Oxzz_temp;
//    A.Oyyy = Oyyy_temp;  A.Oyyz = Oyyz_temp;  A.Oyzz = Oyzz_temp;
//    A.Ozzz = Ozzz_temp;
//}
//
//// Function to subtract multipole `B` from `A` in the reference frame of `A`.
//void subtractMultipole(thread Multipole &A, thread const Multipole &B)
//{
//    // Create a temporary multipole `B_neg` to hold the negated terms.
//    thread Multipole B_neg = B;
//
//    // Negate all components of B.
//    B_neg.M = -B.M;
//
//    B_neg.D = -B.D;
//
//    B_neg.Qxx = -B.Qxx; B_neg.Qxy = -B.Qxy; B_neg.Qxz = -B.Qxz;
//    B_neg.Qyy = -B.Qyy; B_neg.Qyz = -B.Qyz; B_neg.Qzz = -B.Qzz;
//
//    B_neg.Oxxx = -B.Oxxx; B_neg.Oxxy = -B.Oxxy; B_neg.Oxxz = -B.Oxxz;
//    B_neg.Oxyy = -B.Oxyy; B_neg.Oxyz = -B.Oxyz; B_neg.Oxzz = -B.Oxzz;
//    B_neg.Oyyy = -B.Oyyy; B_neg.Oyyz = -B.Oyyz; B_neg.Oyzz = -B.Oyzz;
//    B_neg.Ozzz = -B.Ozzz;
//
//    // Use the existing combineMultipole function to add the negated multipole.
//    combineMultipole(A, B_neg);
//}
//
//
//
//
//
//
//
//// Helper: Q * R. Quadrupole "matrix-vector" multiply.
//inline float3 Q_times_R(thread const Multipole &mp, thread const float3 &R)
//{
//    // Q is:
//    //  [Qxx  Qxy  Qxz]
//    //  [Qxy  Qyy  Qyz]
//    //  [Qxz  Qyz  Qzz]
//    return float3{
//        mp.Qxx*R.x + mp.Qxy*R.y + mp.Qxz*R.z,
//        mp.Qxy*R.x + mp.Qyy*R.y + mp.Qyz*R.z,
//        mp.Qxz*R.x + mp.Qyz*R.y + mp.Qzz*R.z
//    };
//}
//
//// Helper:  T = sum_{i,j,k} O_{ijk} R_i R_j R_k
//inline float octupoleT(thread const Multipole &mp, thread const float3 &R)
//{
//    float x = R.x, y = R.y, z = R.z;
//    return ( mp.Oxxx*x*x*x
//           + mp.Oxxy*x*x*y
//           + mp.Oxxz*x*x*z
//           + mp.Oxyy*x*y*y
//           + mp.Oxyz*x*y*z
//           + mp.Oxzz*x*z*z
//           + mp.Oyyy*y*y*y
//           + mp.Oyyz*y*y*z
//           + mp.Oyzz*y*z*z
//           + mp.Ozzz*z*z*z );
//}
//
//// Helper: grad(T) = partial w.r.t x, y, z of the same expression
//inline float3 gradOctupoleT(thread const Multipole &mp, thread const float3 &R)
//{
//    float x = R.x, y = R.y, z = R.z;
//
//    // dT/dx: differentiate each term that has x in it
//    float dTx = ( 3.f*mp.Oxxx*x*x
//                 + 2.f*mp.Oxxy*x*y
//                 + 2.f*mp.Oxxz*x*z
//                 +    mp.Oxyy*y*y
//                 +    mp.Oxyz*y*z
//                 +    mp.Oxzz*z*z );
//    // dT/dy
//    float dTy = (   mp.Oxxy*x*x
//                 + 2.f*mp.Oxyy*x*y
//                 +    mp.Oxyz*x*z
//                 + 3.f*mp.Oyyy*y*y
//                 + 2.f*mp.Oyyz*y*z
//                 +    mp.Oyzz*z*z );
//    // dT/dz
//    float dTz = (   mp.Oxxz*x*x
//                 +    mp.Oxyz*x*y
//                 + 2.f*mp.Oxzz*x*z
//                 +    mp.Oyyz*y*y
//                 + 2.f*mp.Oyzz*y*z
//                 + 3.f*mp.Ozzz*z*z );
//
//    return float3{ dTx, dTy, dTz };
//}
//
//float3 multipoleAcceleration(thread const Multipole &mp, float3 r, float G)
//{
//    //-------------------------------------
//    // 1) Relative position R = r - R0
//    //-------------------------------------
//    float3 R = r - mp.R0;
//    float R2 = dot(R, R);
//    float R1 = sqrt(R2);
//    if (R1 < 1e-12f) {
//        // Avoid singularities; you might do something safer here
//        return float3{0.f, 0.f, 0.f};
//    }
//    float invR2 = 1.f / R2;
//    float invR  = 1.f / R1;
//    float invR3 = invR2 * invR;  // 1 / R^3
//    float3 accel = float3{0.f, 0.f, 0.f};
//
//    //-------------------------------------
//    // 2) Monopole term: a_M = G*M / R^3 * R
//    //-------------------------------------
//    {
//        float3 aMono = (G * mp.M) * (invR3 * R);
//        accel += aMono;  // remember, you'll do "acc -= accel" in caller if you like
//    }
//
//    //-------------------------------------
//    // 3) Dipole term: a_D = G [ 3(D·R)R / R^5  - D / R^3 ]
//    //-------------------------------------
//    {
//        float dotDR = dot(mp.D, R);
//        float invR5 = invR3 * invR2;  // 1 / R^5
//        float3 aDip = G * ( 3.f * dotDR * R * invR5  -  mp.D * invR3 );
//        accel += aDip;
//    }
//
//    //-------------------------------------
//    // 4) Quadrupole term
//    //
//    //    Q * R => (Qxx*Rx + Qxy*Ry + Qxz*Rz, ...)
//    //    S = R^T Q R = dot(R, Q*R)
//    //
//    //    a_Q = G [ 2(QR)/R^5 - 5 S R / R^7 ]
//    //-------------------------------------
//    {
//        float3 QR = Q_times_R(mp, R);
//        float  S  = dot(R, QR);
//        float invR5 = invR3 * invR2;  // 1 / R^5
//        float invR7 = invR5 * invR2;  // 1 / R^7
//        float3 aQuad = G * ( 2.f * QR * invR5  -  5.f * S * R * invR7 );
//        accel += aQuad;
//    }
//
//    //-------------------------------------
//    // 5) Octupole term
//    //
//    //    T = sum_{i,j,k} O_{ijk} R_i R_j R_k
//    //    ∇T = partial derivatives
//    //
//    //    a_O = G [ ∇T / R^7  -  7 T R / R^9 ]
//    //
//    //-------------------------------------
//    {
//        float Tval = octupoleT(mp, R);
//        float3 gT  = gradOctupoleT(mp, R);
//        float invR7 = invR3 * invR2 * invR2; // 1 / R^7
//        float invR9 = invR7 * invR2;        // 1 / R^9
//
//        float3 aOct = G * ( gT * invR7  -  7.f * Tval * R * invR9 );
//        accel += aOct;
//    }
//
//    //-------------------------------------
//    // Return the total multipole acceleration contribution
//    //-------------------------------------
//    return accel;
//}
//
//void transformMultipole(thread Multipole &mp, thread const float3 &newR0)
//{
//    // Compute the shift vector: Δ = newR0 - oldR0
//    float3 dR = newR0 - mp.R0;
//
//    //---------------------------------------------------------
//    // 1) Monopole: remains the same
//    //---------------------------------------------------------
//    // No changes needed for `mp.M`
//
//    //---------------------------------------------------------
//    // 2) Dipole
//    //---------------------------------------------------------
//    // Shifted dipole: D' = D + M * Δ
//    mp.D += mp.M * dR;
//
//    //---------------------------------------------------------
//    // 3) Quadrupole
//    //---------------------------------------------------------
//    // Shifted quadrupole: Q'_{ij} = Q_{ij} + Δ_i D_j + Δ_j D_i + M Δ_i Δ_j
//    mp.Qxx += dR.x * mp.D.x + dR.x * mp.D.x + mp.M * dR.x * dR.x;
//    mp.Qxy += dR.x * mp.D.y + dR.y * mp.D.x + mp.M * dR.x * dR.y;
//    mp.Qxz += dR.x * mp.D.z + dR.z * mp.D.x + mp.M * dR.x * dR.z;
//    mp.Qyy += dR.y * mp.D.y + dR.y * mp.D.y + mp.M * dR.y * dR.y;
//    mp.Qyz += dR.y * mp.D.z + dR.z * mp.D.y + mp.M * dR.y * dR.z;
//    mp.Qzz += dR.z * mp.D.z + dR.z * mp.D.z + mp.M * dR.z * dR.z;
//
//    //---------------------------------------------------------
//    // 4) Octupole
//    //---------------------------------------------------------
//    // Shifted octupole:
//    // O'_{ijk} = O_{ijk}
//    //            + Δ_i Q_{jk} + Δ_j Q_{ik} + Δ_k Q_{ij}
//    //            + Δ_i Δ_j D_k + Δ_j Δ_k D_i + Δ_k Δ_i D_j
//    //            + M Δ_i Δ_j Δ_k
//
//    // Oxxx
//    mp.Oxxx += dR.x * mp.Qxx + dR.x * mp.Qxx + dR.x * mp.Qxx
//             + dR.x * dR.x * mp.D.x + mp.M * dR.x * dR.x * dR.x;
//
//    // Oxxy
//    mp.Oxxy += dR.x * mp.Qxy + dR.x * mp.Qxy + dR.y * mp.Qxx
//             + dR.x * dR.x * mp.D.y + dR.x * dR.y * mp.D.x
//             + mp.M * dR.x * dR.x * dR.y;
//
//    // Oxxz
//    mp.Oxxz += dR.x * mp.Qxz + dR.x * mp.Qxz + dR.z * mp.Qxx
//             + dR.x * dR.x * mp.D.z + dR.x * dR.z * mp.D.x
//             + mp.M * dR.x * dR.x * dR.z;
//
//    // Oxyy
//    mp.Oxyy += dR.x * mp.Qyy + dR.y * mp.Qxy + dR.y * mp.Qxy
//             + dR.x * dR.y * mp.D.y + dR.y * dR.y * mp.D.x
//             + mp.M * dR.x * dR.y * dR.y;
//
//    // Oxyz
//    mp.Oxyz += dR.x * mp.Qyz + dR.y * mp.Qxz + dR.z * mp.Qxy
//             + dR.x * dR.y * mp.D.z + dR.x * dR.z * mp.D.y + dR.y * dR.z * mp.D.x
//             + mp.M * dR.x * dR.y * dR.z;
//
//    // Oxzz
//    mp.Oxzz += dR.x * mp.Qzz + dR.z * mp.Qxz + dR.z * mp.Qxz
//             + dR.x * dR.z * mp.D.z + dR.z * dR.z * mp.D.x
//             + mp.M * dR.x * dR.z * dR.z;
//
//    // Oyyy
//    mp.Oyyy += dR.y * mp.Qyy + dR.y * mp.Qyy + dR.y * mp.Qyy
//             + dR.y * dR.y * mp.D.y + mp.M * dR.y * dR.y * dR.y;
//
//    // Oyyz
//    mp.Oyyz += dR.y * mp.Qyz + dR.y * mp.Qyz + dR.z * mp.Qyy
//             + dR.y * dR.y * mp.D.z + dR.y * dR.z * mp.D.y
//             + mp.M * dR.y * dR.y * dR.z;
//
//    // Oyzz
//    mp.Oyzz += dR.y * mp.Qzz + dR.z * mp.Qyz + dR.z * mp.Qyz
//             + dR.y * dR.z * mp.D.z + dR.z * dR.z * mp.D.y
//             + mp.M * dR.y * dR.z * dR.z;
//
//    // Ozzz
//    mp.Ozzz += dR.z * mp.Qzz + dR.z * mp.Qzz + dR.z * mp.Qzz
//             + dR.z * dR.z * mp.D.z + mp.M * dR.z * dR.z * dR.z;
//
//    //---------------------------------------------------------
//    // Update the reference origin
//    //---------------------------------------------------------
//    mp.R0 = newR0;
//}
