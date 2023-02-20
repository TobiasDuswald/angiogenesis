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
#include <functional>
#include <iostream>
#include "core/environment/uniform_grid_environment.h"
#include "core/operation/mechanical_forces_op.h"
#include "math.h"
#include "modules/mechanical_forces.h"
#include "modules/tip_cell_finder_operation.h"
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
    vessel_compartment_2->AddBehavior(new SproutingAngiogenesis());
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

double Gaussian(double x, double y, double z) {
  double mu_x = 0.0;
  double mu_y = 0.0;
  double mu_z = 0.0;
  double sigma = 130.0;
  double r = std::sqrt(std::pow(x - mu_x, 2) + std::pow(y - mu_y, 2) +
                       std::pow(z - mu_z, 2));

  return std::exp(-r * r / (2 * sigma * sigma));
}

/// @brief  This function sets up the experiment
/// @param experiment experiment to be set up
/// @param fn Initial concentration of nutrients
/// @param fv Initial concentration of VEGF
/// @param fd Initial concentration of DOX
/// @param ft Initial concentration of TRA
/// @param initialize_random_cells Whether to initialize cells randomly
/// @param initialize_tumor_spheroid Whether to initialize tumor spheroid
/// @param initialize_vasculature Whether to initialize vasculature
/// @param sparam
void SetUpExperiment(const Experiment experiment,
                     std::function<double(double, double, double)>& fn,
                     std::function<double(double, double, double)>& fv,
                     std::function<double(double, double, double)>& fd,
                     std::function<double(double, double, double)>& ft,
                     bool& initialize_random_cells,
                     bool& initialize_tumor_spheroid,
                     bool& initialize_vasculature, const Param* param,
                     const SimParam* sparam) {
  // No use case for random cells at the moment
  initialize_random_cells = false;

  // Computation for Experiment::kVesselsCoupling which is not working as
  // expected in the switch statement below
  const double interval = param->max_bound - param->min_bound;
  const double slope = 1.0 / interval;
  const double offset = -slope * param->min_bound;

  switch (experiment) {
    case Experiment::kAvascularTumorSpheroid:
      initialize_tumor_spheroid = true;
      initialize_vasculature = false;
      break;
    case Experiment::kPorousTumorSpheroid:
      fn = [&sparam](double x, double y, double z) {
        // Compute distance to center (0,0,0)
        double r = std::sqrt(x * x + y * y + z * z);
        if (r > 100) {
          return sparam->initial_concentration_nutrients;
        } else {
          return 0.0;
        }
      };
      initialize_tumor_spheroid = true;
      initialize_vasculature = false;
      break;
    case Experiment::kSpheroidTreatment:
      initialize_tumor_spheroid = true;
      initialize_vasculature = false;
      break;
    case Experiment::kVesselsToCenter:
      fv = Gaussian;
      initialize_tumor_spheroid = false;
      initialize_vasculature = true;
      break;
    case Experiment::kVesselsCoupling:
      fv = [slope, offset](double x, double, double) {
        return slope * x + offset;
      };
      initialize_tumor_spheroid = false;
      initialize_vasculature = true;
      break;
    case Experiment::kSimplifiedGrowth:
      Log::Fatal("SetUpExperiment", "Not implemented yet");
      break;
    case Experiment::kFullScaleModel:
      Log::Fatal("SetUpExperiment", "Not implemented yet");
      break;
    default:
      Log::Fatal("SetUpExperiment", "Unknown experiment");
  };

  // Check if all functions are set
  if (!fn || !fv || !fd || !ft) {
    Log::Fatal("SetUpExperiment", "Not all functions are set");
  }
};

void InitializeVessels(const Experiment experiment, const SimParam* sparam) {
  switch (experiment) {
    case Experiment::kVesselsToCenter:
      PlaceVessel({-200, 0, -400}, {-200, 0, 400},
                  sparam->default_vessel_length);
      PlaceVessel({200, 0, -400}, {200, 0, 400}, sparam->default_vessel_length);
      PlaceVessel({0, -400, 200}, {0, 400, 200}, sparam->default_vessel_length);
      PlaceVessel({0, -400, -200}, {0, 400, -200},
                  sparam->default_vessel_length);
      break;
    case Experiment::kVesselsCoupling:
      PlaceVessel({0, 0, -400}, {0, 0, 400}, sparam->default_vessel_length);
      break;
    case Experiment::kFullScaleModel:
      Log::Fatal("InitializeVessels", "Not implemented yet");
      break;

    default:
      Log::Fatal("InitializeVessels",
                 "No vessel structure defined for this "
                 "experiment");
      break;
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

  constexpr Experiment experiment = Experiment::kVesselsCoupling;

  // ---------------------------------------------------------------------------
  // 1. Define parameters and initialize simulation
  // ---------------------------------------------------------------------------
  auto set_param = [&](Param* param) {
    param->calculate_gradients = true;
    param->visualization_interval =
        param->Get<SimParam>()->visualization_interval /
        param->simulation_time_step;
  };

  // Initialize the simulation
  AngiogenesisSimulation simulation(argc, argv, set_param);
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
  // 2. Get setup for experiment
  // ---------------------------------------------------------------------------
  std::function<double(double, double, double)> initial_nutrient_concentration =
      [&sparam](double, double, double) {
        return sparam->initial_concentration_nutrients;
      };
  std::function<double(double, double, double)> initial_vegf_concentration =
      [&sparam](double, double, double) {
        return sparam->initial_concentration_vegf;
      };
  std::function<double(double, double, double)> initial_dox_concentration =
      [&sparam](double, double, double) {
        return sparam->initial_concentration_dox;
      };
  std::function<double(double, double, double)> initial_tra_concentration =
      [&sparam](double, double, double) {
        return sparam->initial_concentration_tra;
      };
  // Initialize cells randomly at the beginning of the simulation
  bool initialize_random_cells = false;
  // Initialize tumor spheroid at the beginning of the simulation
  bool initialize_tumor_spheroid = false;
  // Initialize vasculature at the beginning of the simulation
  bool initialize_vasculature = false;

  SetUpExperiment(experiment, initial_nutrient_concentration,
                  initial_vegf_concentration, initial_dox_concentration,
                  initial_tra_concentration, initialize_random_cells,
                  initialize_tumor_spheroid, initialize_vasculature, param,
                  sparam);

  // ---------------------------------------------------------------------------
  // 2. Define continuum models for nutrients and VEGF
  // ---------------------------------------------------------------------------

  // Define nutrients with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kNutrients, "Nutrients", sparam->diffusion_nutrients,
      sparam->decay_rate_nutrients, sparam->diffusion_resolution_nutrients);
  ModelInitializer::InitializeSubstance(Substances::kNutrients,
                                        initial_nutrient_concentration);

  // Define VEGF with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kVEGF, "VEGF", sparam->diffusion_vegf,
      sparam->decay_rate_vegf, sparam->diffusion_resolution_vegf);
  ModelInitializer::InitializeSubstance(Substances::kVEGF,
                                        initial_vegf_concentration);

  // Define TRA with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kTRA, "TRA", sparam->diffusion_tra, sparam->decay_rate_tra,
      sparam->diffusion_resolution_tra);
  ModelInitializer::InitializeSubstance(Substances::kTRA,
                                        initial_tra_concentration);

  // Define DOX with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kDOX, "DOX", sparam->diffusion_dox, sparam->decay_rate_dox,
      sparam->diffusion_resolution_dox);
  ModelInitializer::InitializeSubstance(Substances::kDOX,
                                        initial_dox_concentration);

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

    if (initialize_random_cells) {
      // Old single cell initialization
      std::vector<Double3> cell_positions = {
          {0, 50, 0},        {0, 30, 20},      {0, 70, 50},
          {-200, -160, 300}, {-400, -100, 60}, {-300, 100, -200}};
      PlaceTumorCells(cell_positions);
    }

    if (initialize_tumor_spheroid) {
      // Place tumor cells in spheroid
      const uint64_t num_cells = 500;
      const double filled_volume = 0.7;
      const double R =
          std::pow(num_cells * std::pow(sparam->cell_radius, 3) / filled_volume,
                   1.0 / 3.0);
      ModelInitializer::CreateAgentsInSphereRndm({0, 0, 0}, R, num_cells,
                                                 CreateTumorCell);
    }

    if (initialize_vasculature) {
      InitializeVessels(experiment, sparam);
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

  if (sparam->verify_continuum_values) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "VerifyContinuum", OpComputeTarget::kCpu, new VerifyContinuum());
    auto* verify_continuum = NewOperation("VerifyContinuum");
    scheduler->ScheduleOp(verify_continuum, OpType::kPostSchedule);
  }

  // ---------------------------------------------------------------------------
  // 8. Track continuum models
  // ---------------------------------------------------------------------------

  if (sparam->tip_cell_finder_update_frequency <
      std::numeric_limits<int>::max()) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "update tip-cell finder", OpComputeTarget::kCpu,
        new UpdateTipCellFinder());
    auto* update_tip_cell_finder = NewOperation("update tip-cell finder");
    update_tip_cell_finder->frequency_ =
        sparam->tip_cell_finder_update_frequency;
    scheduler->ScheduleOp(update_tip_cell_finder, OpType::kPreSchedule);
  }

  // ---------------------------------------------------------------------------
  // 9. Run simulation and visualize results
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
