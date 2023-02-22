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

#include "random_field.h"
#include "biodynamo.h"

namespace bdm {

RandomField::RandomField(int num_modes, double interval_length,
                         double discretization_goal, double exponent,
                         double max_abs_value, double nu, double sigma,
                         uint32_t random_seed)
    : max_abs_value_(max_abs_value),
      interval_length_(interval_length),
      exponent_(exponent),
      nu_(nu),
      sigma_(sigma),
      num_modes_(num_modes),
      random_seed_(random_seed) {
  // Calculate the number of discretization points
  num_points_ =
      static_cast<int>(std::floor(interval_length / discretization_goal)) + 1;
};

void RandomField::GetRealization(std::vector<double>& random_field) {
  // Resample the random variables
  ResampleRandomVariables();
  // Empty the vector
  random_field.clear();
  // Resize the vector
  random_field.resize(num_points_);
  // Calculate the discretization step
  double discretization_step = interval_length_ / (num_points_ - 1);
  // Track larges absolute value
  double max_abs_value = 0;
  // Calculate the random field
  for (int i = 0; i < num_points_; i++) {
    // Evluate the sinusoidal functions
    double x = i * discretization_step;
    double value = 0;
    for (int j = 0; j < num_modes_; j++) {
      value += amplitudes_[j] *
               std::sin(frequencies_[j] * x / interval_length_ + phases_[j]);
    }
    // Enforce zero beginning and end of the field
    value *= std::pow(std::abs(x * (x - interval_length_)), exponent_);
    // Track the maximum absolute value
    if (std::abs(value) > max_abs_value) {
      max_abs_value = std::abs(value);
    }
    // Add entry to the random field
    random_field[i] = value;
  }
  // Normalize the random field
  for (int i = 0; i < num_points_; i++) {
    random_field[i] *= max_abs_value_ / max_abs_value;
  }
};

void RandomField::ResampleRandomVariables() {
  // Resizing the vectors
  frequencies_.resize(num_modes_);
  amplitudes_.resize(num_modes_);
  phases_.resize(num_modes_);
  // Get random number generator
  Random random;
  random.SetSeed(random_seed_);
  for (int i = 0; i < num_modes_; i++) {
    // Get random frequency
    double frequency = random.Gaus(nu_, sigma_);
    // Get random amplitude
    double amplitude = random.Uniform(-1, 1);
    // Get random phase
    double phase = random.Uniform(0, 2 * Math::kPi);
    // Add to vectors
    frequencies_[i] = frequency;
    amplitudes_[i] = amplitude;
    phases_[i] = phase;
  }
  // Increment the random seed
  random_seed_++;
};

}  // namespace bdm
