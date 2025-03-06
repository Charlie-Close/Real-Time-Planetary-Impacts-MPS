//
//  helpers.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include <metal_stdlib>
#include "poles.h"
using namespace metal;

float integer_powf(const float x, const unsigned int n) {
  switch (n) {
    case 0:
      return 1.f;
    case 1:
      return x;
    case 2:
      return x * x;
    case 3:
      return x * x * x;
    case 4: {
      const float y = x * x;
      return y * y;
    }
    case 5: {
      const float y = x * x;
      return x * y * y;
    }
    case 6: {
      const float y = x * x;
      return y * y * y;
    }
    case 7: {
      const float y = x * x;
      return x * y * y * y;
    }
    case 8: {
      const float y = x * x;
      const float z = y * y;
      return z * z;
    }
    default:
      return pow(x, (float)n);
  }
}

void addPowers(thread Multipole &mp) {
    mp.power[0] = mp.expansion[0];
    mp.power[1] = 0;
    
#if P > 1
    mp.power[2] = mp.expansion[ZZ] * mp.expansion[ZZ];
    mp.power[2] += 0.5 * mp.expansion[YZ] * mp.expansion[YZ];
    mp.power[2] += mp.expansion[YY] * mp.expansion[YY];
    mp.power[2] += 0.5 * mp.expansion[XZ] * mp.expansion[XZ];
    mp.power[2] += 0.5 * mp.expansion[XY] * mp.expansion[XY];
    mp.power[2] += 0.5 * mp.expansion[XX] * mp.expansion[XX];
    mp.power[2] = sqrt(mp.power[2]);
#endif
#if P > 2
    mp.power[3] = mp.expansion[ZZZ] * mp.expansion[ZZZ];
    mp.power[3] += (1 / 3) * mp.expansion[YZZ] * mp.expansion[YZZ];
    mp.power[3] += (1 / 3) * mp.expansion[YYZ] * mp.expansion[YYZ];
    mp.power[3] = mp.expansion[YYY] * mp.expansion[YYY];
    mp.power[3] += (1 / 3) * mp.expansion[XZZ] * mp.expansion[XZZ];
    mp.power[3] += (1 / 6) * mp.expansion[XYZ] * mp.expansion[XYZ];
    mp.power[3] += (1 / 3) * mp.expansion[XYY] * mp.expansion[XYY];
    mp.power[3] += (1 / 3) * mp.expansion[XXZ] * mp.expansion[XXZ];
    mp.power[3] += (1 / 3) * mp.expansion[XXY] * mp.expansion[XXY];
    mp.power[3] = mp.expansion[XXX] * mp.expansion[XXX];
    mp.power[3] = sqrt(mp.power[3]);
#endif
#if P > 3
    mp.power[4] = mp.expansion[ZZZZ] * mp.expansion[ZZZZ];
    mp.power[4] += 2.500000000000000e-01 * mp.expansion[YZZZ] * mp.expansion[YZZZ];
    mp.power[4] += 1.666666666666667e-01 * mp.expansion[YYZZ] * mp.expansion[YYZZ];
    mp.power[4] += 2.500000000000000e-01 * mp.expansion[YYYZ] * mp.expansion[YYYZ];
    mp.power[4] += mp.expansion[YYYY] * mp.expansion[YYYY];
    mp.power[4] += 2.500000000000000e-01 * mp.expansion[XZZZ] * mp.expansion[XZZZ];
    mp.power[4] += 8.333333333333333e-02 * mp.expansion[XYZZ] * mp.expansion[XYZZ];
    mp.power[4] += 8.333333333333333e-02 * mp.expansion[XYYZ] * mp.expansion[XYYZ];
    mp.power[4] += 2.500000000000000e-01 * mp.expansion[XYYY] * mp.expansion[XYYY];
    mp.power[4] += 1.666666666666667e-01 * mp.expansion[XXZZ] * mp.expansion[XXZZ];
    mp.power[4] += 8.333333333333333e-02 * mp.expansion[XXYZ] * mp.expansion[XXYZ];
    mp.power[4] += 1.666666666666667e-01 * mp.expansion[XXYY] * mp.expansion[XXYY];
    mp.power[4] += 2.500000000000000e-01 * mp.expansion[XXXZ] * mp.expansion[XXXZ];
    mp.power[4] += 2.500000000000000e-01 * mp.expansion[XXXY] * mp.expansion[XXXY];
    mp.power[4] += mp.expansion[XXXX] * mp.expansion[XXXX];

    mp.power[4] = sqrt(mp.power[4]);
#endif
}
