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

/// ----------------------------------------------------------------------------
/// This file contains probabilistic state transitions from Rocha et al 2018 /
/// Lima et al 2021. The transitions to dead and proliferative cell states are
/// stochastic and the probabilities depend on certain variables. The functions
/// to compute the probability thresholds are defined in this file.
/// ----------------------------------------------------------------------------

#ifndef TRANSITION_PROBABILITIES_H_
#define TRANSITION_PROBABILITIES_H_

#include "sim_param.h"

namespace bdm {

// Compute the probability for transitioning from quiescent to dead
// (numeric parameter \alpha_{D}(\sigma))
double ComputeProbabilityDeath(const double sigma, const double delta_t,
                               const SimParam* sparam);

// Compute the probability for transitioning from quiescent to
// proliferative (numeric parameter \alpha_{P}(\sigma))
double ComputeProbabilityProliferative(const double sigma, const double delta_t,
                                       const SimParam* sparam);

/// @brief A smooth version of the heaviside step function depending on the
/// value of the concentration. See equation (16) in
/// https://doi.org/10.1016/j.jtbi.2012.02.002. The function is bounded by
/// 1-exp(-alpha * dt) and 1-exp(-(alpha+1) * dt).
/// @param concentration The concentration of the substance
/// @param concentration_threshold The threshold for the substance
/// @param alpha Determines upper and lower bounds of the function.
/// @param k The steepness of transition the function. The higher the value, the
/// steeper the transition. If larger than 0, the function decreases with c, if
/// smaller than 0, the function increases with c.
/// @param delta_t The time step of the simulation (must be small)
/// @return 1-exp(-(alpha+1/(1+exp(2k*(c-c_t))))dt)
inline double SmoothHeavisideForConcentration(double concentration,
                                              double concentration_threshold,
                                              double alpha, double k,
                                              double delta_t) {
  double e = 2 * k * (concentration - concentration_threshold);
  double summand = 1.0 / (1.0 + std::exp(e));
  return 1 - std::exp(-(alpha + summand) * delta_t);
}

}  // namespace bdm

#endif  // TRANSITION_PROBABILITIES_H_
