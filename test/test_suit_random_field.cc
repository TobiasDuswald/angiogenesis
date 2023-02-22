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
#include "util/random_field.h"

namespace bdm {

TEST(RandomField, DiscretizationPoints1) {
  RandomField rf(4, 1, 0.01, 1, 1, 0, 1, 0);
  auto discretization_points = rf.GetNumPoints();
  EXPECT_EQ(discretization_points, 101);
}

TEST(RandomField, DiscretizationPoints2) {
  RandomField rf(4, 10, 1.1, 1, 1, 0, 1, 0);
  auto discretization_points = rf.GetNumPoints();
  EXPECT_EQ(discretization_points, 10);
}

TEST(RandomField, ZeroStartAndEnd) {
  RandomField rf(4, 10, 0.01, 1, 1, 0, 1, 0);
  std::vector<double> random_field;
  rf.GetRealization(random_field);
  ASSERT_EQ(random_field.size(), rf.GetNumPoints());
  EXPECT_DOUBLE_EQ(random_field[0], 0);
  EXPECT_DOUBLE_EQ(random_field[rf.GetNumPoints() - 1], 0);
}

TEST(RandomField, DifferentRealizations) {
  RandomField rf(4, 10, 0.01, 1, 1, 0, 1, 0);
  std::vector<double> random_field1;
  std::vector<double> random_field2;
  rf.GetRealization(random_field1);
  rf.GetRealization(random_field2);
  // Check all points but the first and last one. There is a chance that the
  // we find the same values at intersection points but the probability is
  // very low. We just want to make sure that not all values are the same. E.g.
  // we get different realizations.
  for (int i = 1; i < rf.GetNumPoints() - 1; i++) {
    EXPECT_NE(random_field1[i], random_field2[i]);
  }
}

TEST(RandomField, MaxValue) {
  const double max_abs_value = 1.5;
  RandomField rf(4, 10, 0.01, 1, max_abs_value, 0, 1, 0);
  std::vector<double> random_field;
  rf.GetRealization(random_field);
  double max_abs_value_measured = 0;
  for (int i = 0; i < rf.GetNumPoints(); i++) {
    if (std::abs(random_field[i]) > max_abs_value_measured) {
      max_abs_value_measured = std::abs(random_field[i]);
    }
  }
  EXPECT_DOUBLE_EQ(max_abs_value_measured, max_abs_value);
}

}  // namespace bdm
