//
//  helpers.metal
//  SPH
//
//  Created by Charlie Close on 05/02/2025.
//

#include <metal_stdlib>
#include "poles.h"
using namespace metal;

int binomial(const int n, const int k) {
  /* Hello Pascal! Nice to meet again */
  const int coeffs[9][9] = {
      {1, 0, 0, 0, 0, 0, 0, 0, 0},     {1, 1, 0, 0, 0, 0, 0, 0, 0},
      {1, 2, 1, 0, 0, 0, 0, 0, 0},     {1, 3, 3, 1, 0, 0, 0, 0, 0},
      {1, 4, 6, 4, 1, 0, 0, 0, 0},     {1, 5, 10, 10, 5, 1, 0, 0, 0},
      {1, 6, 15, 20, 15, 6, 1, 0, 0},  {1, 7, 21, 35, 35, 21, 7, 1, 0},
      {1, 8, 28, 56, 70, 56, 28, 8, 1}

  };

  return coeffs[n][k];
}


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
}
