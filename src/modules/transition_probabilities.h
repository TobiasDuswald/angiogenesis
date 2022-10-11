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

// -------------------------------------------------------------------------
// Probabilities for stochastic transitions between cell states
// -------------------------------------------------------------------------

double ComputeProbability_Q_To_SG2(const double nutrients, const double tra,
                                   const double delta_t,
                                   const SimParam* sparam);

double ComputeProbability_Q_To_D(const double nutrients, const double tra,
                                 const double dox, const double delta_t,
                                 const SimParam* sparam);

double ComputeProbability_SG2_To_SG2(const double dox, const double delta_t,
                                     const SimParam* sparam);

double ComputeProbability_SG2_To_D(const double dox, const double delta_t,
                                   const SimParam* sparam);

double ComputeProbability_H_To_D(const double tra, const double dox,
                                 const double delta_t, const SimParam* sparam);

// -------------------------------------------------------------------------
// Individual probabilities
// -------------------------------------------------------------------------

/**** Q -> SG2 ****/

// Compute the probability for transitioning from quiescent Q to proliferative
// SG2 depending on the nutrient concentration
double P_Q_SG2_N(const double nutrients, const double delta_t,
                 const SimParam* sparam);

// Compute the probability for transitioning from quiescent Q to proliferative
// G1 depending on the TRA concentration
double P_Q_SG2_TRA(const double tra, const double delta_t,
                   const SimParam* sparam);

/**** Q -> D ****/

// Compute the probability for transitioning from quiescent Q to proliferative
// SG2 depending on the nutrient concentration
double P_Q_D_N(const double nutrients, const double delta_t,
               const SimParam* sparam);

// Compute the probability for transitioning from quiescent Q to dead D
// depending on the TRA concentration
double P_Q_D_TRA(const double tra, const double delta_t,
                 const SimParam* sparam);

// Compute the probability for transitioning from quiescent Q to dead D
// depending on the DOX concentration
double P_Q_D_DOX(const double dox, const double delta_t,
                 const SimParam* sparam);

/**** SG2 -> G1 ****/

// Compute the probability for reset the proliferative SG2 to the
// dead state depending on the DOX concentration (remains longer in SG2)
double P_SG2_SG2_DOX(const double dox, const double delta_t,
                     const SimParam* sparam);

// Compute the probability for transitioning from proliferative SG2 to the
// dead state depending on the DOX concentration
double P_SG2_D_DOX(const double dox, const double delta_t,
                   const SimParam* sparam);

/**** H -> D ****/

// Compute the probability for transitioning from hypoxic H to dead D depending
// on the DOX concentration
double P_H_D_DOX(const double dox, const double delta_t,
                 const SimParam* sparam);

// Compute the probability for transitioning from hypoxic H to dead D depending
// on the TRA concentration
double P_H_D_TRA(const double tra, const double delta_t,
                 const SimParam* sparam);

// -------------------------------------------------------------------------
// Helper functions
// -------------------------------------------------------------------------

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
/// @param gamma Optional parmater to adjust the function (default = 1)
/// @return 1-exp(-(alpha+gamma/(1+exp(2k*(c-c_t))))dt)
inline double SmoothHeavisideForConcentration(double concentration,
                                              double concentration_threshold,
                                              double alpha, double k,
                                              double delta_t,
                                              double gamma = 1) {
  double e = 2 * k * (concentration - concentration_threshold);
  double summand = gamma / (1.0 + std::exp(e));
  return 1 - std::exp(-(alpha + summand) * delta_t);
}

/// @brief Function to model a linear increase with the concentration.
/// @param concentration The concentration of the substance
/// @param concentration_threshold The threshold for the substance (probability
/// is zero below this value);
/// @param alpha Determines upper bound of the function.
/// @param delta_t The time step of the simulation (must be small)
/// @return 1-exp(-e * dt) with e = max(alpha * (c - c_t)/(1-c_t), 0)
inline double LinearProbabilityIncreaseForConcentration(
    double concentration, double concentration_threshold, double alpha,
    double delta_t) {
  double e = std::max(0.0, alpha * (concentration - concentration_threshold) /
                               (1.0 - concentration_threshold));
  return 1 - std::exp(-e * delta_t);
}

}  // namespace bdm

#endif  // TRANSITION_PROBABILITIES_H_
