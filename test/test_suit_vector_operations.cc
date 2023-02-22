// -----------------------------------------------------------------------------
//
// Copyright (C) 2022 CERN, TUM, and UT Austin. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------

#include <gtest/gtest.h>
#include "util/vector_operations.h"

#define TEST_NAME typeid(*this).name()

namespace bdm {

TEST(VectorOperationsTest, VectorOnUnitCone) {
  // Define parameters
  const double phi = 2.0;
  const double theta = 0.75;
  Double3 z_axis{0, 0, 1};

  // Get vector in unit sphere
  Double3 my_vector = VectorOnUnitCone(phi, theta);

  // Comute references
  const double z = std::cos(theta);
  const double r_squared = 1 - std::pow(z, 2);

  // Checks
  EXPECT_DOUBLE_EQ(z, my_vector[2]);
  EXPECT_DOUBLE_EQ(1.0, my_vector.Norm());
  EXPECT_DOUBLE_EQ(r_squared,
                   std::pow(my_vector[0], 2) + std::pow(my_vector[1], 2));
  EXPECT_DOUBLE_EQ(theta, std::acos(z_axis * my_vector));
  EXPECT_DOUBLE_EQ(phi - Math::kPi, std::atan(my_vector[1] / my_vector[0]));
}

TEST(VectorOperationsTest, VectorOnConeAroundAxis) {
  // Define parameters
  const double phi = 2.0;
  const double theta = 0.75;
  Double3 rot_axis{1.23, -0.75, 2.0};
  rot_axis.Normalize();

  // Get vector in unit sphere
  Double3 my_vector = VectorOnConeAroundAxis(rot_axis, phi, theta);

  // Checks if the parameter theta is realized correctly, phi is ignored because
  // it is usually rotationally invariant in our use cases.
  EXPECT_DOUBLE_EQ(theta, std::acos(rot_axis * my_vector));
  EXPECT_DOUBLE_EQ(1.0, my_vector.Norm());
}

TEST(VectorOperationsTest, VectorOnConeAroundZAxis) {
  // Define parameters
  const double phi = 2.0;
  const double theta = 0.75;
  Double3 rot_axis{0, 0, 1};

  // Get vector in unit sphere
  Double3 my_vector = VectorOnConeAroundAxis(rot_axis, phi, theta);

  // Checks if the parameter theta is realized correctly
  double phi_computed = std::atan(my_vector[1] / my_vector[0]);
  if (my_vector[0] < 0) {
    phi_computed += Math::kPi;
  }
  EXPECT_DOUBLE_EQ(theta, std::acos(rot_axis * my_vector));
  EXPECT_DOUBLE_EQ(phi, phi_computed);
  EXPECT_DOUBLE_EQ(1.0, my_vector.Norm());

  rot_axis = {0, 0, -1};
  my_vector = VectorOnConeAroundAxis(rot_axis, phi, theta);
  phi_computed = std::atan(my_vector[1] / my_vector[0]);
  if (my_vector[0] < 0) {
    phi_computed += Math::kPi;
  }
  EXPECT_DOUBLE_EQ(theta, std::acos(rot_axis * my_vector));
  EXPECT_DOUBLE_EQ(1.0, my_vector.Norm());
}

TEST(VectorOperationsTest, GetOrthogonalSystem) {
  std::vector<Double3> vectors{{1, 0, 0}, {0, 1, 0},  {0, 0, 1},
                               {1, 1, 1}, {1, 1, 0},  {1, 0, 1},
                               {0, 1, 1}, {1, 1, -1}, {1.23, -0.75, 2.0}};

  // Push back 10 random vectors
  Random random;
  for (int i = 0; i < 10; i++) {
    vectors.push_back(random.UniformArray<3>(-1, 1));
  }

  for (const auto& a : vectors) {
    // Define parameters
    // Double3 a{1.23, -0.75, 2.0};
    // a.Normalize();
    Double3 b{0, 0, 0};
    Double3 c{0, 0, 0};

    // Get orthogonal system
    GetOrthogonalSystem(a, b, c);

    // Checks
    EXPECT_DOUBLE_EQ(1.0, b.Norm());
    EXPECT_DOUBLE_EQ(1.0, c.Norm());
    EXPECT_LT(std::abs(a * b), 1e-14);
    EXPECT_LT(std::abs(a * c), 1e-14);
    EXPECT_LT(std::abs(b * c), 1e-14);
  }
}

}  // namespace bdm
