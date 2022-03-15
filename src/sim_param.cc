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

#include "sim_param.h"
#include <cmath>

namespace bdm {

double SimParam::ComputeProbabilityDeath(const double sigma,
                                         const double delta_t) const {
  double intensity{apoptosis_rate};
  intensity += gamma * 1 / (1 + std::exp(2 * k * (sigma - hypoxic_threshold)));
  return 1 - std::exp(-intensity * delta_t);
}

double SimParam::ComputeProbabilityProliferative(const double sigma,
                                                 const double delta_t) const {
  double intensity = std::max(qp_transition_rate * (sigma - hypoxic_threshold) /
                                  (1 - hypoxic_threshold),
                              0.0);
  return 1 - std::exp(-intensity * delta_t);
}

}  // namespace bdm
