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

#include "treatment.h"
#include <cmath>

namespace bdm {

bool Treatment::IsTraApplied(double t) const {
  if (t >= tra_start_1_ && t < tra_end_1_) {
    return true;
  }
  if (t >= tra_start_2_ && t < tra_end_2_) {
    return true;
  }
  return false;
}

bool Treatment::IsDoxApplied(double t) const {
  if (t >= dox_start_ && t < dox_end_) {
    return true;
  }
  return false;
};

/// ODE for the vessel permeability.
double Treatment::VesselPermeabilityODE(double x, double t) const {
  if (IsTraApplied(t)) {
    return (max_vessel_permeability_ - x) / vessel_permeability_growth_;
  } else {
    return -x / vessel_permeability_decay_;
  }
}

/// Precompute the vessel permeability.
void Treatment::PrecomputeVesselPermeability(double t_end, double time_step,
                                             double time_step_ode) {
  // Initial condition
  vessel_permeability.clear();
  double x = vessel_permeability_0_;

  // Determine the number of backup steps, i.e. the number of steps that are
  // written to the vessel permeability vector.
  int backup_steps = std::ceil(t_end / time_step);

  // Adapt ode time step to backup steps
  int num_ode_steps = std::ceil(time_step / time_step_ode);
  time_step_ode = time_step / num_ode_steps;

  // Precompute vessel permeability
  double time = 0;
  vessel_permeability.push_back(x);
  for (int i = 0; i < backup_steps; i++) {
    for (int j = 0; j < num_ode_steps; j++) {
      // Forward Euler
      double tmp = VesselPermeabilityODE(x, time);
      x += time_step_ode * tmp;
      // Avoid overshooting
      if (x < 0) {
        x = 0;
      }
      if (x > max_vessel_permeability_) {
        x = max_vessel_permeability_;
      }
      // Update time
      time += time_step_ode;
    }
    vessel_permeability.push_back(x);
  }
};

}  // namespace bdm
