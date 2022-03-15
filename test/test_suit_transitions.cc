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
#include "modules/transition_probabilities.h"
#include "sim_param.h"

namespace bdm {

// The SimParam struct can compute the probability of cell to transition into
// the death state. This computation in tested here.
TEST(SimParam, ComputeProbabilityDeath) {
  auto sparam = SimParam();
  // Set values used in reference solution.
  sparam.apoptosis_rate = 0.000408 / 60.0;
  sparam.k = 50.0;
  sparam.gamma = 0.0245 / 60.0;
  sparam.hypoxic_threshold = 0.0538;

  // Test some numeric values of the function
  EXPECT_FLOAT_EQ(4.100812398e-06,
                  ComputeProbabilityDeath(0.01, 0.01, &sparam));
  EXPECT_FLOAT_EQ(3.8054246159e-06,
                  ComputeProbabilityDeath(0.03, 0.01, &sparam));
  EXPECT_FLOAT_EQ(2.4929787294e-06,
                  ComputeProbabilityDeath(0.05, 0.01, &sparam));
  EXPECT_FLOAT_EQ(7.425862775e-07,
                  ComputeProbabilityDeath(0.07, 0.01, &sparam));
  EXPECT_FLOAT_EQ(1.0783971216e-07,
                  ComputeProbabilityDeath(0.10, 0.01, &sparam));
  EXPECT_FLOAT_EQ(0.00816809991, ComputeProbabilityDeath(0.01, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.00758197442, ComputeProbabilityDeath(0.03, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.00497355438, ComputeProbabilityDeath(0.05, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.00148407078, ComputeProbabilityDeath(0.07, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.00021565617, ComputeProbabilityDeath(0.10, 20.0, &sparam));
}

// The SimParam struct can compute the probability of cell to transition into
// the proliferative state. This computation in tested here.
TEST(SimParam, ComputeProbabilityProliferative) {
  auto sparam = SimParam();
  // Set values used in reference solution.
  sparam.qp_transition_rate = 0.0493 / 60.0;
  sparam.hypoxic_threshold = 0.0538;

  // Test some numeric values of the function
  EXPECT_FLOAT_EQ(0.0, ComputeProbabilityProliferative(0.01, 0.01, &sparam));
  EXPECT_FLOAT_EQ(0.0, ComputeProbabilityProliferative(0.03, 0.01, &sparam));
  EXPECT_FLOAT_EQ(0.0, ComputeProbabilityProliferative(0.05, 0.01, &sparam));
  EXPECT_FLOAT_EQ(1.4067849363197382e-07,
                  ComputeProbabilityProliferative(0.07, 0.01, &sparam));
  EXPECT_FLOAT_EQ(4.011941702186661e-07,
                  ComputeProbabilityProliferative(0.10, 0.01, &sparam));
  EXPECT_FLOAT_EQ(0.0, ComputeProbabilityProliferative(0.01, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.0, ComputeProbabilityProliferative(0.03, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.005994663318752647,
                  ComputeProbabilityProliferative(0.4, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.01629904273146321,
                  ComputeProbabilityProliferative(1.0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.049881690944160284,
                  ComputeProbabilityProliferative(3.0, 20.0, &sparam));
  EXPECT_FLOAT_EQ(0.09811829568804986,
                  ComputeProbabilityProliferative(6.0, 20.0, &sparam));
}
}  // namespace bdm
