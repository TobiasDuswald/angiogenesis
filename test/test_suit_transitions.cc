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
#include <array>
#include "modules/transition_probabilities.h"
#include "sim_param.h"

namespace bdm {

TEST(SimParam, ComputeProbability_Q_To_D) {
  auto sparam = SimParam();
  // Set values used in reference solution.
  sparam.alpha_Q_D_N = 0.000408 / 60.0;
  sparam.k_Q_D_N = 50.0;
  sparam.gamma_Q_D_N = 0.0245 / 60.0;
  sparam.threshold_Q_D_N = 0.0538;

  // Test some numeric values of the function
  EXPECT_FLOAT_EQ(4.100812398e-06,
                  ComputeProbability_Q_To_D(0.01, 0, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(3.8054246159e-06,
                  ComputeProbability_Q_To_D(0.03, 0, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(2.4929787294e-06,
                  ComputeProbability_Q_To_D(0.05, 0, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(7.425862775e-07,
                  ComputeProbability_Q_To_D(0.07, 0, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(1.0783971216e-07,
                  ComputeProbability_Q_To_D(0.10, 0, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(0.00816809991,
                  ComputeProbability_Q_To_D(0.01, 0, 0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.00758197442,
                  ComputeProbability_Q_To_D(0.03, 0, 0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.00497355438,
                  ComputeProbability_Q_To_D(0.05, 0, 0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.00148407078,
                  ComputeProbability_Q_To_D(0.07, 0, 0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.00021565617,
                  ComputeProbability_Q_To_D(0.10, 0, 0, 20.0, &sparam));
}

TEST(SimParam, ComputeProbability_Q_To_SG2) {
  auto sparam = SimParam();
  // Set values used in reference solution.
  sparam.alpha_Q_SG2_N = 0.0493 / 60.0;
  sparam.threshold_Q_SG2_N = 0.0538;

  // Test some numeric values of the function
  EXPECT_FLOAT_EQ(0.0, ComputeProbability_Q_To_SG2(0.01, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(0.0, ComputeProbability_Q_To_SG2(0.03, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(0.0, ComputeProbability_Q_To_SG2(0.05, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(1.4067849363197382e-07,
                  ComputeProbability_Q_To_SG2(0.07, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(4.011941702186661e-07,
                  ComputeProbability_Q_To_SG2(0.10, 0, 0.01, &sparam));
  EXPECT_FLOAT_EQ(0.0, ComputeProbability_Q_To_SG2(0.01, 0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.0, ComputeProbability_Q_To_SG2(0.03, 0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.005994663318752647,
                  ComputeProbability_Q_To_SG2(0.4, 0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.01629904273146321,
                  ComputeProbability_Q_To_SG2(1.0, 0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.049881690944160284,
                  ComputeProbability_Q_To_SG2(3.0, 0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.09811829568804986,
                  ComputeProbability_Q_To_SG2(6.0, 0, 20.0, &sparam));
}

TEST(TransitionHelperFunctions, SmoothHeavisideForConcentration) {
  // Parameter values used in reference solution
  const double k = 8.0;
  const double alpha = 3.0;
  const double c_t = 0.5;
  const double dt = 0.1;

  // Concentrations (x values) used for testing
  std::array<double, 13> concentrations(
      {0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0});

  // Reference values for the sigmoid function via Wolfram Alpha
  std::array<double, 13> expected_results(
      {0.329657, 0.329569, 0.329133, 0.327049, 0.318325, 0.308573, 0.295312,
       0.281797, 0.271522, 0.262078, 0.259786, 0.259305, 0.259207});

  // Test the sigmoid function for all concentrations
  for (size_t i = 0; i < concentrations.size(); i++) {
    // Reference solution is only given with 6 digits precision, thus we demand
    // the same precision for the test -> 1e-6
    EXPECT_LT(std::abs(expected_results[i] -
                       SmoothHeavisideForConcentration(concentrations[i], c_t,
                                                       alpha, k, dt)),
              1e-6);
  }
}

TEST(TransitionHelperFunctions, LinearProbabilityIncreaseForConcentration) {
  // Parameter values used in reference solution
  const double alpha = 3.0;
  const double c_t = 0.5;
  const double dt = 0.1;

  // Concentrations (x values) used for testing
  std::array<double, 13> concentrations(
      {0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0});

  // Reference values for the sigmoid function via Python
  std::array<double, 13> expected_results({0., 0., 0., 0., 0., 0., 0.,
                                           0.02955447, 0.05823547, 0.11307956,
                                           0.16472979, 0.21337214, 0.25918178});

  // Test the function for all concentrations
  for (size_t i = 0; i < concentrations.size(); i++) {
    // Reference solution is only given with 6 digits precision, thus we demand
    // the same precision for the test -> 1e-6
    EXPECT_LT(std::abs(expected_results[i] -
                       LinearProbabilityIncreaseForConcentration(
                           concentrations[i], c_t, alpha, dt)),
              1e-6);
  }
}

}  // namespace bdm
