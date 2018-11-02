// This file implements functions to calculate various physical quantities
#ifndef CS_UTIL_H
#define CS_UTIL_H

#include "common.h"
#include <cstddef>
#include <valarray>
#include <vector>

namespace cs {

// Calculate total energy
FLOAT CALC_ENERGY(const valarray<FLOAT> &mass,
                  const vector<valarray<FLOAT>> &pos,
                  const vector<valarray<FLOAT>> &vel)
{
  FLOAT k = 0.0;  // kinetic energy
  for (size_t i = 0; i < mass.size(); i++) k += 0.5 * mass[i] * (vel[i] * vel[i]).sum();
  FLOAT u = 0.0; // potential energy
  for (size_t i = 0; i < mass.size() - 1; i++) {
    for (size_t j = i + 1; j < mass.size(); j++) {
      auto r = pos[i] - pos[j];
      u -= (G * mass[i] * mass[j]) / sqrt((r * r).sum());
    }
  }
  return k + u;
}


// Calculate total angular momentum
valarray<FLOAT> CALC_ANGULAR_MOMENTUM(const valarray<FLOAT> &mass,
                                      const vector<valarray<FLOAT>> &pos,
                                      const vector<valarray<FLOAT>> &vel)
{
  valarray<FLOAT> h(DIM);
  for (size_t i = 0; i < mass.size(); i++) {
    h += mass[i] * (pos[i].cshift(1) * vel[i].cshift(-1) - pos[i].cshift(-1) * vel[i].cshift(1));
  }
  return h;
}


}  // namespace cs

#endif //CS_UTIL_H
