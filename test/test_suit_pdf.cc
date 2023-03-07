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
#include "biodynamo.h"
#include "core/util/random.h"
#include "util/pdf.h"

#define TEST_NAME typeid(*this).name()

namespace bdm {

inline void CompareVectors(const std::vector<double>& v1,
                           const std::vector<double>& v2, double tol) {
  ASSERT_EQ(v1.size(), v2.size());
  for (int i = 0; i < v1.size(); i++) {
    EXPECT_NEAR(v1[i], v2[i], tol);
  }
}

TEST(PDF, GEV) {
  // Define x parameters to test
  const std::vector<double> x_values = {0.0,  1.0,  2.0,  3.0, 4.0,  5.0,
                                        6.0,  7.0,  8.0,  9.0, 10.0, 11.0,
                                        12.0, 13.0, 14.0, 15.0};

  // --------------------------------------------------------------------------
  // Test Basic Functionality
  // --------------------------------------------------------------------------

  // Define reference values
  double mu = 0.0;
  double sigma = 1.0;
  double xi = 0;
  std::vector<double> y_val = {
      3.67879441e-01, 2.54646380e-01, 1.18204952e-01, 4.73690097e-02,
      1.79832297e-02, 6.69269968e-03, 2.47261557e-03, 9.11050816e-04,
      3.35350112e-04, 1.23394575e-04, 4.53978687e-05, 1.67014218e-05,
      6.14417460e-06, 2.26032430e-06, 8.31528028e-07, 3.05902227e-07};

  // Calculate values with function
  std::vector<double> y_values;
  for (auto x : x_values) {
    y_values.push_back(gev_pdf(x, mu, sigma, xi));
  }

  // Compare the values
  CompareVectors(y_values, y_val, 1e-8);

  // --------------------------------------------------------------------------
  // Test xi = -0.5
  // --------------------------------------------------------------------------

  // Define reference values
  mu = 0.0;
  sigma = 1.0;
  xi = -0.5;
  y_val = {0.36787944, 0.18997937, 0.0973501,  0.0545372,
           0.0331422,  0.02149529, 0.01467833, 0.01044518,
           0.00768632, 0.00581507, 0.0045028,  0.00355616,
           0.00285656, 0.0023286,  0.00192284, 0.00160595};

  // Calculate values with function
  y_values.clear();
  for (auto x : x_values) {
    y_values.push_back(gev_pdf(x, mu, sigma, xi));
  }

  // // Compare the values
  CompareVectors(y_values, y_val, 1e-8);

  // --------------------------------------------------------------------------
  // Test location shift
  // --------------------------------------------------------------------------

  // Define reference values
  mu = -3.0;
  sigma = 1.0;
  xi = -0.5;

  y_val = {0.0545372,  0.0331422,  0.02149529, 0.01467833,
           0.01044518, 0.00768632, 0.00581507, 0.0045028,
           0.00355616, 0.00285656, 0.0023286,  0.00192284,
           0.00160595, 0.00135491, 0.0011535,  0.00099005};

  // Calculate values with function
  y_values.clear();
  for (auto x : x_values) {
    y_values.push_back(gev_pdf(x, mu, sigma, xi));
  }

  // Compare the values
  CompareVectors(y_values, y_val, 1e-8);

  // --------------------------------------------------------------------------
  // Test scale parameter
  // --------------------------------------------------------------------------

  // Define reference values
  mu = 0.0;
  sigma = 2.0;
  xi = -0.5;

  y_val = {0.18393972, 0.13498686, 0.09498969, 0.0673047,
           0.04867505, 0.03602765, 0.0272686,  0.02106418,
           0.0165711,  0.01324962, 0.01074764, 0.00883066,
           0.00733916, 0.00616253, 0.00522259, 0.00446314};

  // Calculate values with function
  y_values.clear();
  for (auto x : x_values) {
    y_values.push_back(gev_pdf(x, mu, sigma, xi));
  }

  // Compare the values
  CompareVectors(y_values, y_val, 1e-8);
}

TEST(PDF, WALD) {
  // Define x parameters to test
  const std::vector<double> x_values = {0.0,  1.0,  2.0,  3.0, 4.0,  5.0,
                                        6.0,  7.0,  8.0,  9.0, 10.0, 11.0,
                                        12.0, 13.0, 14.0, 15.0};

  // --------------------------------------------------------------------------
  // Test Basic Functionality
  // --------------------------------------------------------------------------

  // Define reference values
  double mu = 0.0;
  double sigma = 1.0;
  std::vector<double> y_val = {
      0.00000000e+00, 3.98942280e-01, 1.09847822e-01, 3.94183580e-02,
      1.61896995e-02, 7.20416893e-03, 3.37989353e-03, 1.64628783e-03,
      8.24609311e-04, 4.22073556e-04, 2.19794800e-04, 1.16079415e-04,
      6.20254871e-05, 3.34712362e-05, 1.82154793e-05, 9.98579307e-06};

  // Calculate values with function
  std::vector<double> y_values;
  for (auto x : x_values) {
    y_values.push_back(wald_pdf(x, mu, sigma));
  }

  // Compare the values
  CompareVectors(y_values, y_val, 1e-8);

  // --------------------------------------------------------------------------
  // Test location shift
  // --------------------------------------------------------------------------

  // Define reference values
  mu = -3.0;
  sigma = 1.0;

  y_val = {3.94183580e-02, 1.61896995e-02, 7.20416893e-03, 3.37989353e-03,
           1.64628783e-03, 8.24609311e-04, 4.22073556e-04, 2.19794800e-04,
           1.16079415e-04, 6.20254871e-05, 3.34712362e-05, 1.82154793e-05,
           9.98579307e-06, 5.50930754e-06, 3.05671327e-06, 1.70444000e-06};

  // Calculate values with function
  y_values.clear();
  for (auto x : x_values) {
    y_values.push_back(wald_pdf(x, mu, sigma));
  }

  // Compare the values
  CompareVectors(y_values, y_val, 1e-8);

  // --------------------------------------------------------------------------
  // Test scale parameter
  // --------------------------------------------------------------------------

  // Define reference values
  mu = 0.0;
  sigma = 2.0;

  y_val = {0.,         0.43939129, 0.19947114, 0.09989689,
           0.05492391, 0.03217641, 0.01970918, 0.01247427,
           0.00809485, 0.0053572,  0.00360208, 0.0024538,
           0.00168995, 0.00117474, 0.00082314, 0.0005808};

  // Calculate values with function
  y_values.clear();
  for (auto x : x_values) {
    y_values.push_back(wald_pdf(x, mu, sigma));
  }

  // Compare the values
  CompareVectors(y_values, y_val, 1e-8);
}

TEST(PDF, GEV_RNG) {
  // Define GEV parameters
  double mu = 7.53401234807571;
  double sigma = 2.3153691046974907;
  double xi = -0.3292535751714517;
  // Expected mean and standard deviation [via Scipy]
  // >>> import scipy.stats as st
  // >>> st.genextreme.mean(-0.3292535751714517,
  //                        7.53401234807571,
  //                        2.3153691046974907)
  // >>> st.genextreme.std(-0.3292535751714517,
  //                        7.53401234807571,
  //                        2.3153691046974907)
  double expected_mean = 9.973398037878336;
  double expected_std_dev = 6.282287338421184;

  // Define a user defined random number generator for GEV
  auto distribution = [](const double* x, const double* params) {
    return gev_pdf(x[0], params[0], params[1], params[2]);
  };

  // Random number generator
  Simulation simulation("test");
  auto* random = simulation.GetRandom();
  auto rng =
      random->GetUserDefinedDistRng1D(distribution, {mu, sigma, xi}, 0, 500);

  // Generate 10000 random numbers
  std::vector<double> random_numbers;
  for (int i = 0; i < 1e7; i++) {
    random_numbers.push_back(rng.Sample());
  }

  // Calculate the mean and standard deviation
  double mean = 0.0;
  double std_dev = 0.0;
  for (auto x : random_numbers) {
    mean += x;
  }
  mean /= random_numbers.size();
  for (auto x : random_numbers) {
    std_dev += (x - mean) * (x - mean);
  }
  std_dev = std::sqrt(std_dev / random_numbers.size());

  // Compare the values
  EXPECT_NEAR(mean, expected_mean, 0.1);
  EXPECT_NEAR(std_dev, expected_std_dev, 0.15);
}

TEST(PDF, WALD_RNG) {
  // Define GEV parameters
  double mu = 10.977589465183868;
  double sigma = 63.81303941790211;
  // Expected mean and standard deviation [via Scipy]
  // >>> import scipy.stats as st
  // >>> st.wald.mean(10.977589465183868, 63.81303941790211)
  // >>> st.wald.std(10.977589465183868, 63.81303941790211)
  double expected_mean = 74.79062888308599;
  double expected_std_dev = 63.81303941790211;

  // Define a user defined random number generator for GEV
  auto distribution = [](const double* x, const double* params) {
    return wald_pdf(x[0], params[0], params[1]);
  };

  // Random number generator
  Simulation simulation("test");
  auto* random = simulation.GetRandom();
  auto rng =
      random->GetUserDefinedDistRng1D(distribution, {mu, sigma}, 0, 5000);

  // Generate 10000 random numbers
  std::vector<double> random_numbers;
  for (int i = 0; i < 1e7; i++) {
    random_numbers.push_back(rng.Sample());
  }

  // Calculate the mean and standard deviation
  double mean = 0.0;
  double std_dev = 0.0;
  for (auto x : random_numbers) {
    mean += x;
  }
  mean /= random_numbers.size();
  for (auto x : random_numbers) {
    std_dev += (x - mean) * (x - mean);
  }
  std_dev = std::sqrt(std_dev / random_numbers.size());

  // Compare the values
  EXPECT_NEAR(mean, expected_mean, 1.0);
  EXPECT_NEAR(std_dev, expected_std_dev, 1.5);
}

}  // namespace bdm
