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

#include "angiogenesis_simulation.h"
#include <fstream>
#include <iostream>
#include "core/environment/uniform_grid_environment.h"
#include "core/operation/mechanical_forces_op.h"
#include "math.h"
#include "modules/mechanical_forces.h"
#include "modules/tumor_cell.h"
#include "sim_param.h"
#include "util/timeseries_counters.h"
#include "util/visualize.h"

namespace bdm {

// Initialize parameter group Uid, part of the BioDynaMo API, needs to be part
// of a cc file, depends on #include "sim-param.h". With this, we can access the
// simulation parameters anywhere in the simulation.
const ParamGroupUid SimParam::kUid = ParamGroupUidGenerator::Get()->NewUid();

auto CreateTumorCell(const Double3& position) {
  // Connect to active simulation, get parameters and a random generator
  auto* sim = Simulation::GetActive();
  auto* random = sim->GetRandom();
  auto* param = sim->GetParam();
  auto* sparam = param->Get<SimParam>();
  // Create cells at a random position (randomness through model initializer)
  // that can divide and are quiescent
  int cell_state = CellState::kHypoxic;
  auto* tumor_cell = new TumorCell(position, cell_state);
  // Set (random) radius, nuclear radius, and action radius
  double random_radius =
      random->Gaus(sparam->cell_radius, sparam->cell_radius_sigma);
  tumor_cell->SetActionRadiusFactor(sparam->action_radius_factor);
  tumor_cell->SetRadii(random_radius, sparam->cell_nuclear_radius,
                       sparam->action_radius_factor * random_radius);
  // Cells gain half their volume during the growth phase.
  double growth_rate = 2.0 / 3.0 * Math::kPi * pow(random_radius, 3) /
                       (sparam->duration_growth_phase);
  tumor_cell->SetGrowthRate(growth_rate);
  // Add a Secretion behaviour to the TumorCell such that it can secrete VEGF.
  tumor_cell->AddBehavior(new Secretion(
      "VEGF", sparam->secretion_rate_vegf * param->simulation_time_step));
  return tumor_cell;
}

void PlaceTumorCells() {
  auto* tumor_cell1 = CreateTumorCell({0, 5, 0});
  auto* tumor_cell2 = CreateTumorCell({0, 3, 2});
  auto* tumor_cell3 = CreateTumorCell({0, 7, 5});
  auto* rm = Simulation::GetActive()->GetResourceManager();
  rm->AddAgent(tumor_cell1);
  rm->AddAgent(tumor_cell2);
  rm->AddAgent(tumor_cell3);
}

// -----------------------------------------------------------------------------
// MAIN SIMULATION
// -----------------------------------------------------------------------------
int Simulate(int argc, const char** argv) {
  // Register the simulation parameter
  Param::RegisterParamGroup(new SimParam());

  // ---------------------------------------------------------------------------
  // 1. Define parameters and initialize simulation
  // ---------------------------------------------------------------------------
  auto set_param = [&](Param* param) {
    param->statistics = true;
    param->calculate_gradients = false;
    param->visualization_interval =
        param->Get<SimParam>()->visualization_interval /
        param->simulation_time_step;
  };

  // Initialize the simulation
  Simulation simulation(argc, argv, set_param);
  // Get a pointer to the resource manager
  const auto* rm = simulation.GetResourceManager();
  // Get a pointer to the param object
  const auto* param = simulation.GetParam();
  // Get a pointer to an instance of SimParam
  const auto* sparam = param->Get<SimParam>();
  // Get a pointer to the scheduler
  auto* scheduler = simulation.GetScheduler();
  // Get a pointer to the environment
  auto* env =
      dynamic_cast<UniformGridEnvironment*>(simulation.GetEnvironment());

  // ---------------------------------------------------------------------------
  // 2. Define continuum models for nutrients and VEGF
  // ---------------------------------------------------------------------------

  // Define nutrients with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kNutrients, "Nutrients", sparam->diffusion_nutrients,
      sparam->decay_rate_nutrients, sparam->diffusion_resolution);
  auto SetInitialValuesGridNutrients = [&](double x, double y, double z) {
    return sparam->initial_nutrient_concentration;
  };
  ModelInitializer::InitializeSubstance(Substances::kNutrients,
                                        SetInitialValuesGridNutrients);

  // Define nutrients with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kVEGF, "VEGF", sparam->diffusion_vegf,
      sparam->decay_rate_vegf, sparam->diffusion_resolution);
  auto SetInitialValuesGridVEGF = [&](double x, double y, double z) {
    return 0.0;
  };
  ModelInitializer::InitializeSubstance(Substances::kNutrients,
                                        SetInitialValuesGridVEGF);

  // ---------------------------------------------------------------------------
  // 3. Define initial configurations of agents
  // ---------------------------------------------------------------------------

  PlaceTumorCells();

  // ---------------------------------------------------------------------------
  // 4. Track simulation results over time with timeseries objects
  // ---------------------------------------------------------------------------

  // Collect the number of Cells in different states over time
  auto* ts = simulation.GetTimeSeries();
  ts->AddCollector("quiescent_cells", CountQuiescent, GetSimulatedTime);
  ts->AddCollector("G1_cells", CountG1, GetSimulatedTime);
  ts->AddCollector("SG2_cells", CountSG2, GetSimulatedTime);
  ts->AddCollector("hypoxic_cells", CountHypoxic, GetSimulatedTime);
  ts->AddCollector("dead_cells", CountDead, GetSimulatedTime);

  // ---------------------------------------------------------------------------
  // 5. Use force module typically used by UT Austin
  // ---------------------------------------------------------------------------

  // Use custom force module implemented in MechanicalInteractionForce
  auto* custom_force = new MechanicalInteractionForce(
      sparam->adhesion_scale_parameter, sparam->repulsive_scale_parameter);
  auto* op = scheduler->GetOps("mechanical forces")[0];
  auto* force_implementation = op->GetImplementation<MechanicalForcesOp>();
  force_implementation->SetInteractionForce(custom_force);

  // ---------------------------------------------------------------------------
  // 6. Specific fix for force and environment combination
  // ---------------------------------------------------------------------------

  // Set box length manually because our interaction range is larger than the
  // cell's diameter.
  // Todo(tobias): check if factor 2 in necessary and how it influences the
  // computational runtime
  env->SetBoxLength(static_cast<int32_t>(
      std::ceil(2 * sparam->action_radius_factor *
                (sparam->cell_radius + 5 * sparam->cell_radius_sigma))));

  // ---------------------------------------------------------------------------
  // 7. Run simulation and visualize results
  // ---------------------------------------------------------------------------

  // Finalize initialization
  scheduler->FinalizeInitialization();

  // Test if correct number of Agents were initialized
  std::cout << "Agents im Simulation: " << rm->GetNumAgents() << "\n";

  // Run simulation for a defined number of timesteps
  u_int64_t time_steps{static_cast<u_int64_t>(
      ceil(sparam->total_sim_time / param->simulation_time_step))};
  scheduler->Simulate(time_steps);
  std::cout << "Simulation completed successfully!" << std::endl;

  // Save timeseries data
  PlotAndSaveTimeseries();

  return 0;
}

}  // namespace bdm