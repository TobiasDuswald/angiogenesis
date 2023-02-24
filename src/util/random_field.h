// --------------------------------------------------------------------------
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
// --------------------------------------------------------------------------

#ifndef RANDOM_FIELD_H_
#define RANDOM_FIELD_H_

#include <cmath>
#include <vector>

namespace bdm {

/// This class generates a random field. In the current version, we simple use
/// random sinusoidal functions to generate the field.
class RandomField {
 public:
  RandomField() = default;

  /// @brief Constructor
  /// @param num_modes Number of sinusoidal functions
  /// @param interval_length Length of the interval [a, b] -> (b - a)
  /// @param discretization_goal Anticipated discretization step (h); may
  ///                            deviate from this value
  /// @param exponent Exponent used to smooth the beginning and end of the
  ///                 random field (x * (x - interval_length))^exponent
  /// @param max_abs_value Largest absolute value of the random field used for
  ///                      normalization
  /// @param nu Mean of the frequency distribution
  /// @param sigma Standard deviation of the frequency distribution
  /// @param random_seed Random seed
  RandomField(int num_modes, double interval_length, double discretization_goal,
              double exponent, double max_abs_value, double nu, double sigma,
              unsigned int random_seed);

  /// Get a random field realization
  void GetRealization(std::vector<double>& random_field);

  /// Get the number of discretization points
  int GetNumPoints() const { return num_points_; }

 private:
  void ResampleRandomVariables();

  // The random vector for the frequencies
  std::vector<double> frequencies_{};
  // The random vector for the amplitudes
  std::vector<double> amplitudes_{};
  // The random vector for the phases
  std::vector<double> phases_{};
  // Maximum value of the random field
  double max_abs_value_ = 0;
  // Interval length of the random field
  double interval_length_ = 0;
  // Exponent on enforcing of the zero beginning and end of the field
  double exponent_ = 0;
  // Mean of the frequency distribution
  double nu_ = 0;
  // Standard deviation of the frequency distribution
  double sigma_ = 0;
  // Discretization points of the random field
  int num_points_ = 0;
  // The number of modes
  int num_modes_ = 0;
  // Random seed
  unsigned int random_seed_ = 0;
};

}  // namespace bdm

#endif  // RANDOM_FIELD_H_
