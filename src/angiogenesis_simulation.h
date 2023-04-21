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

#ifndef ANGIOGENESIS_SIMULATION_H_
#define ANGIOGENESIS_SIMULATION_H_

#include "biodynamo.h"
#include "modules/tipcell_finder.h"
#include "modules/treatment.h"

namespace bdm {

// This function contains the core simulation code. It creates the agents in
// an environment and simulates the system for multiple timesteps.
int Simulate(int argc, const char** argv);

class AngiogenesisSimulation : public Simulation {
 public:
  AngiogenesisSimulation(int argc, const char** argv,
                         const std::function<void(Param*)>& set_param)
      : Simulation(argc, argv, set_param){};

  explicit AngiogenesisSimulation(const std::string& simulation_name)
      : Simulation(simulation_name){};

  const TipCellFinder* GetTipCellFinder() const { return &tip_cell_finder_; }
  void UpdateTipCellFinder() { tip_cell_finder_.Update(); }

  const std::vector<double>& GetVesselPermeability() const {
    return treatment_.GetVesselPermeability();
  }

  double GetVesselPermeability(size_t time_step) const {
    return treatment_.GetVesselPermeability()[time_step];
  }

  void PrecomputeVesselPermeability(double t_end, double time_step,
                                    double time_step_ode) {
    treatment_.PrecomputeVesselPermeability(t_end, time_step, time_step_ode);
  }

  Treatment& GetTreatment() { return treatment_; }

 private:
  TipCellFinder tip_cell_finder_;
  Treatment treatment_;
};

}  // namespace bdm

#endif  // ANGIOGENESIS_SIMULATION_H_
