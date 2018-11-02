// This file implements test for "converter.h"
#include "../converter.h"
#include "gtest/gtest.h"
#include <valarray>
#include <cmath>

namespace cs {


TEST(ConvertOelemToCartesianTest, EccentricityEqualZero)
{
  valarray<double> mass {1.0, 1.0e-7};
  cs::OrbitalElement oelem(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  vector<cs::OrbitalElement> oelems {oelem};
  vector<valarray<double>> pos(2, valarray<double>(3));
  vector<valarray<double>> vel(2, valarray<double>(3));
  cs::ConvertOelemToCartesian(mass, oelems, &pos, &vel);
  for (auto row : pos) for (auto x : row) EXPECT_TRUE(isnan(x));
  for (auto row : vel) for (auto x : row) EXPECT_TRUE(isnan(x));
}


TEST(InterconversionTest, FiveBody)
{
  const size_t num = 5;  // sun + planet
  valarray<double> mass {1.0, 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7};
  cs::OrbitalElement oelem0(1.0, 1.0E-2, 0.0, M_PI/2.0, -M_PI/2.0, 0.0);
  cs::OrbitalElement oelem1(2.0, 2.0E-2, 1.0, 2.0, -1.0, 1.0);
  cs::OrbitalElement oelem2(3.0, 5.0E-2, 1.5, 1.0, 2.0, 2.0);
  cs::OrbitalElement oelem3(4.0, 1.0E-1, 2.0, 0.5, 0.5, 3.0);
  vector<cs::OrbitalElement> init_oelems {oelem0, oelem1, oelem2, oelem3};
  vector<valarray<double>> pos(num, valarray<double>(DIM));
  vector<valarray<double>> vel(num, valarray<double>(DIM));
  cs::ConvertOelemToCartesian(mass, init_oelems, &pos, &vel);
  vector<cs::OrbitalElement> ret_oelems = cs::ConvertCartesianToOelem(mass, pos, vel);
  for (size_t i = 0; i < init_oelems.size(); i++) {
    EXPECT_NEAR(init_oelems[i].sa, ret_oelems[i].sa, 1.0e-3);
    EXPECT_NEAR(init_oelems[i].oe, ret_oelems[i].oe, 1.0e-3);
    EXPECT_NEAR(fmod(init_oelems[i].peri - ret_oelems[i].peri, 2.0 * M_PI), 0.0, 1.0e-3);
    EXPECT_NEAR(init_oelems[i].incl, ret_oelems[i].incl, 1.0e-3);
    EXPECT_NEAR(fmod(init_oelems[i].node - ret_oelems[i].node, 2.0 * M_PI), 0.0, 1.0e-3);
    EXPECT_NEAR(fmod(init_oelems[i].l - ret_oelems[i].l, 2.0 * M_PI), 0.0, 1.0e-3);
  }
}

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
