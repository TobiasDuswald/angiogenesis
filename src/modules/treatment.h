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

#ifndef TREATMENT_H_
#define TREATMENT_H_

#include <vector>

namespace bdm {

class Treatment {
 public:
  Treatment() = default;
  ~Treatment() = default;

  /// Indicator function for the TRA treatment. Returns true if the treatment is
  /// applied, false otherwise.
  /// @param t - current simulation time [min]
  bool IsTraApplied(double t) const;

  /// Indicator function for the DOX treatment. Returns true if the treatment is
  /// applied, false otherwise.
  /// @param t - current simulation time [min]
  bool IsDoxApplied(double t) const;

  /// ODE for the vessel permeability.
  double VesselPermeabilityODE(double x, double t) const;

  /// Precompute the vessel permeability.
  void PrecomputeVesselPermeability(double t_end, double time_step,
                                    double time_step_ode);

  const std::vector<double>& GetVesselPermeability() const {
    return vessel_permeability;
  }

  // Set the six treatment parameters.
  void SetTreatmentParameters(double tra_start_1, double tra_end_1,
                              double tra_start_2, double tra_end_2,
                              double dox_start, double dox_end) {
    tra_start_1_ = tra_start_1;
    tra_end_1_ = tra_end_1;
    tra_start_2_ = tra_start_2;
    tra_end_2_ = tra_end_2;
    dox_start_ = dox_start;
    dox_end_ = dox_end;
  }

 private:
  std::vector<double> vessel_permeability;

  // Treatment parameters. Start end in minutes.
  double tra_start_1_ = 102 * 60 * 24;
  double tra_end_1_ = 103 * 60 * 24;
  double tra_start_2_ = 105 * 60 * 24;
  double tra_end_2_ = 106 * 60 * 24;
  double dox_start_ = 108 * 60 * 24;
  double dox_end_ = 109 * 60 * 24;

  // Vessel permeability parameters
  double vessel_permeability_0_ = 0.0;
  double max_vessel_permeability_ = 1;
  double vessel_permeability_decay_ = 10.0 * 60 * 24;
  double vessel_permeability_growth_ = 0.4 * 60 * 24;
};

}  // namespace bdm

#endif  // TREATMENT_H_
