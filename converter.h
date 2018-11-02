// This file implements functions to convert between orbital elements(sa, oe, peri, incl, node, l) and cartesian coordinates(x, y, z, vx, vy, vz).
//
// - sa: semimajor axis
// - oe: orbital eccentricity
// - peri: argument of periapsis
// - incl: inclination
// - node: longitude of the ascending node
// - l: mean anomaly
#ifndef CS_CONVERTER_H
#define CS_CONVERTER_H

#include "common.h"
#include "oelem.h"
#include <cmath>
#include <numeric>
#include <valarray>
#include <vector>

namespace cs {

// Inline function
template <class T>
inline T SQR(T x) {return x * x;}
template <class T>
inline T CUBE(T x) {return x * x * x;}

// Convert cartesian coordinates to orbital elements
// Range:
//   incl: [0, PI]
//   peri: [0, 2*PI]
//   l: [0, 2*PI]
vector<OrbitalElement> ConvertCartesianToOelem(const valarray<FLOAT> &mass,
                                               const vector<valarray<FLOAT>> &pos,
                                               const vector<valarray<FLOAT>> &vel)
{
  vector<OrbitalElement> oelems(pos.size() - 1);
  for (size_t i = 0; i < pos.size() - 1; i++) {
    OrbitalElement tmp_oelem;

    // Relative position and relative velocity
    valarray<FLOAT> rel_pos = pos[i+1] - pos[0];
    valarray<FLOAT> rel_vel = vel[i+1] - vel[0];
    FLOAT rel_pos_norm = sqrt((rel_pos * rel_pos).sum());
    FLOAT rel_vel_norm = sqrt((rel_vel * rel_vel).sum());

    FLOAT mu = G * (mass[0] + mass[i+1]);

    // Semi-major axis
    tmp_oelem.sa = mu / (2.0 * mu / rel_pos_norm - rel_vel_norm * rel_vel_norm);

    // Angular momentum
    valarray<FLOAT> h = rel_pos.cshift(1) * rel_vel.cshift(-1) - rel_pos.cshift(-1) * rel_vel.cshift(1);
    FLOAT h_norm = sqrt((h*h).sum());

    // Inclination [0, PI]
    tmp_oelem.incl = acos((rel_pos[0] * rel_vel[1] - rel_pos[1] * rel_vel[0]) / h_norm);

    // Longitude of the ascending node
    tmp_oelem.node = atan2(rel_pos[1] * rel_vel[2] - rel_pos[2] * rel_vel[1],
                           rel_pos[0] * rel_vel[2] - rel_pos[2] * rel_vel[0]);

    // Position and velocity
    // (XY plane match with orbital plane and the direction of the X axis match with the direction of ascending node)
    valarray<FLOAT> tmp_pos {rel_pos[0] * cos(tmp_oelem.node) + rel_pos[1] * sin(tmp_oelem.node),
                             (-rel_pos[0] * sin(tmp_oelem.node) + rel_pos[1] * cos(tmp_oelem.node))
                             * cos(tmp_oelem.incl) + rel_pos[2] * sin(tmp_oelem.incl)};
    valarray<FLOAT> tmp_vel {rel_vel[0] * cos(tmp_oelem.node) + rel_vel[1] * sin(tmp_oelem.node),
                             (-rel_vel[0] * sin(tmp_oelem.node) + rel_vel[1] * cos(tmp_oelem.node))
                             * cos(tmp_oelem.incl) + rel_vel[2] * sin(tmp_oelem.incl)};

    // Eccentric vector
    valarray<FLOAT> e_vec {h_norm  * tmp_vel[1] / mu - tmp_pos[0] / rel_pos_norm,
                           -h_norm * tmp_vel[0] / mu - tmp_pos[1] / rel_pos_norm};

    // Eccentricity
    tmp_oelem.oe = sqrt((e_vec*e_vec).sum());

    // Argument of periapsis [0, 2*PI]
    tmp_oelem.peri = atan2(e_vec[1], e_vec[0]);
    if (tmp_oelem.peri < 0.0) tmp_oelem.peri += 2.0 * M_PI;

    // True anomaly [0, 2*PI]
    FLOAT f = atan2(tmp_pos[1], tmp_pos[0]) - tmp_oelem.peri;
    f = fmod(f, 2.0 * M_PI);
    if (f < 0.0) f += 2.0 * M_PI;

    // Eccentric anomaly [0, 2*PI]
    FLOAT u = atan2((rel_pos_norm * sin(f)) / (tmp_oelem.sa * sqrt(1.0 - SQR(tmp_oelem.oe))),
                    (rel_pos_norm * cos(f)) / tmp_oelem.sa + tmp_oelem.oe );
    if (u < 0.0) u += 2.0 * M_PI;

    // Mean anomaly [0, 2*PI]
    tmp_oelem.l = u - tmp_oelem.oe * sin(u);
    if (tmp_oelem.l < 0.0) tmp_oelem.l += 2.0 * M_PI;

    // Update oelems
    oelems[i] = tmp_oelem;
  }
  return oelems;
}



// Convert orbital elements to cartesian coordinates
void ConvertOelemToCartesian(const valarray<FLOAT> &mass,
                             const vector<OrbitalElement> &oelems,
                             vector<valarray<FLOAT>> *p_pos,
                             vector<valarray<FLOAT>> *p_vel)
{
  // Exception handling
  for (OrbitalElement oelem : oelems) {
    if (oelem.oe == 0.0) {
      for (size_t i = 0; i < oelems.size() + 1; i++) {
        (*p_pos)[i] = numeric_limits<FLOAT>::quiet_NaN();
        (*p_vel)[i] = numeric_limits<FLOAT>::quiet_NaN();
      }
      return;
    }
  }

  (*p_pos)[0] = 0.0;
  (*p_vel)[0] = 0.0;
  for (size_t i = 0; i < oelems.size(); i++) {
    FLOAT mu = G * (mass[0] + mass[i+1]);
    // Mean motion
    FLOAT n = sqrt(mu / CUBE(oelems[i].sa));
    // Eccentric anomaly
    FLOAT u = oelems[i].l;
    FLOAT u_last;
    // Permissible error
    constexpr FLOAT MARGIN = 1.0e-20;
    for (;;) {
      u_last = u;
      u -= (u - oelems[i].oe * sin(u) - oelems[i].l) / (1.0 - oelems[i].oe * cos(u));
      if (abs(u - u_last) < MARGIN) break;
    }
    // Calculate position and velocity under the following conditions.
    // 1. orbital plane agrees with xy plane
    // 2. suppose peri = 0
    (*p_pos)[i+1] = {oelems[i].sa * (cos(u) - oelems[i].oe),
                     oelems[i].sa * sqrt(1.0 - oelems[i].oe * oelems[i].oe) * sin(u),
                     0.0};
    (*p_vel)[i+1] = {(-oelems[i].sa * n * sin(u)) / (1.0 - oelems[i].oe * cos(u)),
                     oelems[i].sa * n * sqrt(1.0 - oelems[i].oe * oelems[i].oe) * cos(u)
                     / (1.0 - oelems[i].oe * cos(u)),
                     0.0};
    // Rotate position and velocity (-peri) around the Z-axis
    (*p_pos)[i+1] = {(*p_pos)[i+1][0] * cos(oelems[i].peri) - (*p_pos)[i+1][1] * sin(oelems[i].peri),
                     (*p_pos)[i+1][0] * sin(oelems[i].peri) + (*p_pos)[i+1][1] * cos(oelems[i].peri),
                     0.0};
    (*p_vel)[i+1]= {(*p_vel)[i+1][0] * cos(oelems[i].peri) - (*p_vel)[i+1][1] * sin(oelems[i].peri),
                    (*p_vel)[i+1][0] * sin(oelems[i].peri) + (*p_vel)[i+1][1] * cos(oelems[i].peri),
                    0.0};
    if (oelems[i].incl != 0.0) {
      // Rotate position and velocity (-incl) around the X-axis, and then rotate them (-node) around the Z-axis.
      (*p_pos)[i+1] = {(*p_pos)[i+1][0] * cos(oelems[i].node)
                       - (*p_pos)[i+1][1] * sin(oelems[i].node) * cos(oelems[i].incl)
                       + (*p_pos)[i+1][2] * sin(oelems[i].node) * sin(oelems[i].incl),
                       (*p_pos)[i+1][0] * sin(oelems[i].node)
                       + (*p_pos)[i+1][1] * cos(oelems[i].node) * cos(oelems[i].incl)
                       - (*p_pos)[i+1][2] * cos(oelems[i].node) * sin(oelems[i].incl),
                       (*p_pos)[i+1][1] * sin(oelems[i].incl) + (*p_pos)[i+1][2] * cos(oelems[i].incl)};
      (*p_vel)[i+1] = {(*p_vel)[i+1][0] * cos(oelems[i].node)
                       - (*p_vel)[i+1][1] * sin(oelems[i].node) * cos(oelems[i].incl)
                       + (*p_vel)[i+1][2] * sin(oelems[i].node) * sin(oelems[i].incl),
                       (*p_vel)[i+1][0] * sin(oelems[i].node)
                       + (*p_vel)[i+1][1] * cos(oelems[i].node) * cos(oelems[i].incl)
                       - (*p_vel)[i+1][2] * cos(oelems[i].node) * sin(oelems[i].incl),
                       (*p_vel)[i+1][1] * sin(oelems[i].incl) + (*p_vel)[i+1][2] * cos(oelems[i].incl)};
    }
  }
  // Adjust barycenter position and barycenter velocity
  valarray<FLOAT> pos_moment = {0.0, 0.0, 0.0};
  valarray<FLOAT> vel_moment = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < oelems.size() + 1; i++) {
    pos_moment += mass[i] * (*p_pos)[i];
    vel_moment += mass[i] * (*p_vel)[i];
  }
  valarray<FLOAT> bc_pos = pos_moment / mass.sum();
  valarray<FLOAT> bc_vel = vel_moment / mass.sum();
  for (size_t i = 0; i < oelems.size() + 1; i++) {
    (*p_pos)[i] -= bc_pos;
    (*p_vel)[i] -= bc_vel;
  }
}


}  // namespace cs


#endif // CS_CONVERTER_H