// This file implements functions to write cartesian coordinates data (position and velocity) and orbital elements data in CSV format.
#ifndef CS_CSVWRITER_H
#define CS_CSVWRITER_H

#include "common.h"
#include "oelem.h"
#include <fstream>
#include <limits>

namespace cs {

// Write CSV header of cartesian coordinates data
void WritePosVelHeader(size_t num, ofstream *p_ofs) {
  *p_ofs << "t";
  for (size_t i = 0; i < num; i++) *p_ofs << ",x" << i << ",y" << i << ",z" << i;
  for (size_t i = 0; i < num; i++) *p_ofs << ",vx" << i << ",vy" << i << ",vz" << i;
  *p_ofs << endl;
}


// Write data of cartesian coordinates
template <class F, class CTR2D_F>
void WritePosVelData(FLOAT t, const vector<valarray<FLOAT>> &pos, const vector<valarray<FLOAT>> &vel, ofstream *p_ofs)
{
  constexpr int d = numeric_limits<float>::max_digits10;
  *p_ofs << setprecision(d) << scientific << t;
  for (valarray<FLOAT> row : pos) for (FLOAT x : row) *p_ofs << "," << x;
  for (valarray<FLOAT> row : vel) for (FLOAT x : row) *p_ofs << "," << x;
  *p_ofs << endl;
}


// Write CSV header of orbital elements data
void WriteOelemHeader(size_t num, ofstream *p_ofs) {
  *p_ofs << "t";
  for (size_t i = 0; i < num; i++) {
    *p_ofs << ",sa" << i << ",oe" << i << ",peri" << i << ",incl" << i << ",node" << i << ",l" << i;
  }
  *p_ofs << endl;
}


// Write data of orbital elements
void WriteOelemData(FLOAT t, const vector<OrbitalElement> &oelems, ofstream *p_ofs)
{
  constexpr int d = numeric_limits<float>::max_digits10;
  *p_ofs << setprecision(d) << scientific << t;
  for (auto oelem : oelems) {
    *p_ofs << "," << oelem.sa << "," << oelem.oe << "," << oelem.peri << "," << oelem.incl << "," << oelem.node
           << "," << oelem.l;
  }
  *p_ofs << endl;
}


}  // namespace cs

#endif  // CS_CSVWRITER_H
