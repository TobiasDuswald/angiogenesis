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
#include "modules/vessel.h"
#include "neuroscience/neuroscience.h"
#include "sim_param.h"
#include "util/analysis.h"
#include "util/timeseries_counters.h"

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
  // Add the continuum interactions to the tumor cell.
  tumor_cell->AddBehavior(new PointContinuumInteraction(
      sparam->nutrient_consumption_rate_tcell, sparam->vegf_supply_rate_tcell,
      sparam->dox_consumption_rate_tcell, sparam->tra_consumption_rate_tcell));
  // Add cell cycle
  tumor_cell->AddBehavior(new ProgressInCellCycle());
  return tumor_cell;
}

void PlaceTumorCells(std::vector<Double3>& positions) {
  auto* rm = Simulation::GetActive()->GetResourceManager();
  for (auto pos : positions) {
    auto* tumor_cell = CreateTumorCell(pos);
    rm->AddAgent(tumor_cell);
  }
}

void inline PlaceVessel(Double3 start, Double3 end, double compartment_length) {
  auto* rm = Simulation::GetActive()->GetResourceManager();
  auto* param = Simulation::GetActive()->GetParam();
  auto* sparam = param->Get<SimParam>();

  // Compute parameters for straight line between start and end.
  Double3 direction = end - start;
  double distance = direction.Norm();
  direction.Normalize();
  int n_compartments = std::floor(distance / compartment_length);

  // Warn if chosen parameters are not selected ideally
  if (abs(n_compartments * compartment_length - distance) > 1e-2) {
    Log::Warning("PlaceVessel", "Vessel will be shorter than expected.");
  }

  // The setup requires us to define a NeuronSoma, which is kind of a left over
  // from the neuroscience module.
  const Double3 tmp = start;
  auto* soma = new neuroscience::NeuronSoma(tmp);
  rm->AddAgent(soma);

  // Define a first neurite
  Vessel v;  // Used for prototype argument (virtual+template not supported c++)
  auto* vessel_compartment_1 =
      dynamic_cast<Vessel*>(soma->ExtendNewNeurite(direction, &v));
  vessel_compartment_1->SetPosition(start +
                                    direction * compartment_length * 0.5);
  vessel_compartment_1->SetMassLocation(start + direction * compartment_length);
  vessel_compartment_1->SetActualLength(compartment_length);
  vessel_compartment_1->SetDiameter(15);
  vessel_compartment_1->ProhibitGrowth();

  Vessel* vessel_compartment_2{nullptr};
  for (int i = 1; i < n_compartments; i++) {
    // Compute location of next vessel element
    Double3 agent_position =
        start + direction * compartment_length * (static_cast<double>(i) + 0.5);
    Double3 agent_end_position =
        start + direction * compartment_length * (static_cast<double>(i) + 1.0);
    // Create new vessel
    vessel_compartment_2 = new Vessel();
    // Set position an length
    vessel_compartment_2->SetPosition(agent_position);
    vessel_compartment_2->SetMassLocation(agent_end_position);
    vessel_compartment_2->SetActualLength(compartment_length);
    vessel_compartment_2->SetRestingLength(compartment_length);
    vessel_compartment_2->SetSpringAxis(direction);
    vessel_compartment_2->SetDiameter(15);
    vessel_compartment_2->ProhibitGrowth();
    // Add behaviours
    // vessel_compartment_2->AddBehavior(new SproutingAngiogenesis());
    vessel_compartment_2->AddBehavior(new ApicalGrowth());
    vessel_compartment_2->AddBehavior(new LineContinuumInteraction(
        sparam->nutrient_supply_rate_vessel,
        sparam->vegf_consumption_rate_vessel, sparam->dox_supply_rate_vessel,
        sparam->tra_supply_rate_vessel));
    // Add Agent to the resource manager
    rm->AddAgent(vessel_compartment_2);
    // Connect vessels (AgentPtr API is currently bounded to base
    // classes but this is a 'cosmetic' problem)
    vessel_compartment_1->SetDaughterLeft(
        vessel_compartment_2->GetAgentPtr<neuroscience::NeuriteElement>());
    vessel_compartment_2->SetMother(
        vessel_compartment_1->GetAgentPtr<neuroscience::NeuronOrNeurite>());
    std::swap(vessel_compartment_1, vessel_compartment_2);
  }
}

// -----------------------------------------------------------------------------
// MAIN SIMULATION
// -----------------------------------------------------------------------------
int Simulate(int argc, const char** argv) {
  // Register the simulation parameter
  Param::RegisterParamGroup(new SimParam());
  neuroscience::InitModule();

  // ---------------------------------------------------------------------------
  // 1. Define parameters and initialize simulation
  // ---------------------------------------------------------------------------
  auto set_param = [&](Param* param) {
    param->calculate_gradients = true;
    param->visualization_interval =
        param->Get<SimParam>()->visualization_interval /
        param->simulation_time_step;
    param->visualization_compress_pv_files = true;
  };

  // Initialize the simulation
  Simulation simulation(argc, argv, set_param);
  // Get a pointer to the resource manager
  auto* rm = simulation.GetResourceManager();
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
      sparam->decay_rate_nutrients, sparam->diffusion_resolution_nutrients);
  auto SetInitialValuesGridNutrients = [&sparam](double x, double y, double z) {
    // return sparam->initial_concentration_nutrients;
    // return sparam->hypoxic_threshold * 0.5;
    return sparam->hypoxic_threshold * 2;
  };
  ModelInitializer::InitializeSubstance(Substances::kNutrients,
                                        SetInitialValuesGridNutrients);

  // Define VEGF with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kVEGF, "VEGF", sparam->diffusion_vegf,
      sparam->decay_rate_vegf, sparam->diffusion_resolution_vegf);
  auto SetInitialValuesGridVEGF = [&sparam](double x, double y, double z) {
    return sparam->initial_concentration_vegf;
  };
  ModelInitializer::InitializeSubstance(Substances::kVEGF,
                                        SetInitialValuesGridVEGF);

  // Define TRA with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kTRA, "TRA", sparam->diffusion_tra, sparam->decay_rate_tra,
      sparam->diffusion_resolution_tra);
  auto SetInitialValuesGridTRA = [&sparam](double x, double y, double z) {
    return sparam->initial_concentration_tra;
  };
  ModelInitializer::InitializeSubstance(Substances::kTRA,
                                        SetInitialValuesGridTRA);

  // Define DOX with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kDOX, "DOX", sparam->diffusion_dox, sparam->decay_rate_dox,
      sparam->diffusion_resolution_dox);
  auto SetInitialValuesGridDOX = [&sparam](double x, double y, double z) {
    return sparam->initial_concentration_dox;
  };
  ModelInitializer::InitializeSubstance(Substances::kDOX,
                                        SetInitialValuesGridDOX);

  // Define upper and lower threshold for nutrients
  rm->ForEachDiffusionGrid([&](DiffusionGrid* grid) {
    grid->SetUpperThreshold(1.0);
    grid->SetLowerThreshold(0.0);
  });

  // ---------------------------------------------------------------------------
  // 3. Define initial configurations of agents
  // ---------------------------------------------------------------------------
  {
    Timing timer_set_up("Initialize agents");

    if (sparam->initialize_random_cells) {
      // Old single cell initialization
      std::vector<Double3> cell_positions = {
          {0, 50, 0},        {0, 30, 20},      {0, 70, 50},
          {-200, -160, 300}, {-400, -100, 60}, {-300, 100, -200}};
      PlaceTumorCells(cell_positions);
    }

    if (sparam->initialize_tumor_spheroid) {
      // Place tumor cells in spheroid
      const uint64_t num_cells = 500;
      const double filled_volume = 0.7;
      const double R =
          std::pow(num_cells * std::pow(sparam->cell_radius, 3) / filled_volume,
                   1.0 / 3.0);
      ModelInitializer::CreateAgentsInSphereRndm({0, 0, 0}, R, num_cells,
                                                 CreateTumorCell);
    }

    if (sparam->initialize_vasculature) {
      // Place vessels
      PlaceVessel({-200, 0, -400}, {-200, 0, 400},
                  sparam->default_vessel_length);
      PlaceVessel({200, 0, -400}, {200, 0, 400}, sparam->default_vessel_length);
      PlaceVessel({0, -400, 200}, {0, 400, 200}, sparam->default_vessel_length);
      PlaceVessel({0, -400, -200}, {0, 400, -200},
                  sparam->default_vessel_length);
    }
  }

  // ---------------------------------------------------------------------------
  // 4. Track simulation results over time with timeseries objects
  // ---------------------------------------------------------------------------

  // Collect the number of Cells in different states over time
  DefineAndRegisterCollectors();

  // ---------------------------------------------------------------------------
  // 5. Use force module typically used by UT Austin
  // ---------------------------------------------------------------------------

  // Use custom force module implemented in MechanicalInteractionForce
  // Note that the force module currently does not support any forces
  // between vessels and cells.
  auto* custom_force = new MechanicalInteractionForce(
      sparam->adhesion_scale_parameter, sparam->repulsive_scale_parameter);
  auto* op = scheduler->GetOps("mechanical forces")[0];
  auto* force_implementation = op->GetImplementation<MechanicalForcesOp>();
  force_implementation->SetInteractionForce(custom_force);

  // ---------------------------------------------------------------------------
  // 6. Specific fix for force and environment combination
  // ---------------------------------------------------------------------------

  // Set box length manually because our interaction range is larger than the
  // cell's diameter. In the current setup we restrict vessel growth once we
  // come close to a tumor cell. We set the UniformGrid to larger box sizes
  // such that we can resolve more long distance relationships than with the
  // UniformGrid.
  double distance_for_growth_stop = 60;
  double box_length =
      std::ceil(2 * sparam->action_radius_factor *
                (sparam->cell_radius + 5 * sparam->cell_radius_sigma));
  box_length = std::max(box_length, distance_for_growth_stop);
  env->SetBoxLength(static_cast<int32_t>(box_length));

  // ---------------------------------------------------------------------------
  // 7. Track continuum models
  // ---------------------------------------------------------------------------

  OperationRegistry::GetInstance()->AddOperationImpl(
      "VerifyContinuum", OpComputeTarget::kCpu, new VerifyContinuum());
  auto* verify_continuum = NewOperation("VerifyContinuum");
  scheduler->ScheduleOp(verify_continuum, OpType::kPostSchedule);

  // ---------------------------------------------------------------------------
  // 8 . Run simulation and visualize results
  // ---------------------------------------------------------------------------

  // Finalize initialization
  scheduler->FinalizeInitialization();
  scheduler->PrintInfo(std::cout);

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
