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

}  // namespace bdm

#endif  // TRANSITION_PROBABILITIES_H_
