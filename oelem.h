// This file implements "Orbital Elements" struct.
#ifndef CS_OELEM_H
#define CS_OELEM_H

#include "common.h"
#include <limits>

namespace cs {

// Orbital element
struct OrbitalElement {
  OrbitalElement()
          : sa(numeric_limits<FLOAT>::quiet_NaN()),
            oe(numeric_limits<FLOAT>::quiet_NaN()),
            peri(numeric_limits<FLOAT>::quiet_NaN()),
            incl(numeric_limits<FLOAT>::quiet_NaN()),
            node(numeric_limits<FLOAT>::quiet_NaN()),
            l(numeric_limits<FLOAT>::quiet_NaN()) {}

  OrbitalElement(FLOAT sa, FLOAT oe, FLOAT peri, FLOAT incl, FLOAT node, FLOAT l)
          : sa(sa),
            oe(oe),
            peri(peri),
            incl(incl),
            node(node),
            l(l) {}

  FLOAT sa;    // semimajor axis
  FLOAT oe;    // orbital eccentricity
  FLOAT peri;  // argument of periapsis
  FLOAT incl;  // inclination
  FLOAT node;  // longitude of the ascending node
  FLOAT l;     // mean anomaly
};

}


#endif //CS_OELEM_H
