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

#include "transition_probabilities.h"
#include <cmath>

namespace bdm {

double ComputeProbabilityDeath(const double sigma, const double delta_t,
                               const SimParam* sparam) {
  return SmoothHeavisideForConcentration(sigma, sparam->hypoxic_threshold,
                                         sparam->apoptosis_rate, sparam->k,
                                         delta_t, sparam->gamma);
}

double ComputeProbabilityProliferative(const double sigma, const double delta_t,
                                       const SimParam* sparam) {
  double intensity = std::max(sparam->qp_transition_rate *
                                  (sigma - sparam->hypoxic_threshold) /
                                  (1 - sparam->hypoxic_threshold),
                              0.0);
  return 1 - std::exp(-intensity * delta_t);
}

}  // namespace bdm
