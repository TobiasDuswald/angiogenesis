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
#include <Math/Integrator.h>
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
#include "util/data_parser.h"
#include "util/pdf.h"
#include "util/random_field.h"
#include "util/timeseries_counters.h"
#include "util/vector_operations.h"

namespace bdm {

// Initialize parameter group Uid, part of the BioDynaMo API, needs to be part
// of a cc file, depends on #include "sim-param.h". With this, we can access the
// simulation parameters anywhere in the simulation.
const ParamGroupUid SimParam::kUid = ParamGroupUidGenerator::Get()->NewUid();

// -----------------------------------------------------------------------------
// Helper function (structure; implementation below main simulation)
// -----------------------------------------------------------------------------

/// Function to create a TumorCell with random properties. This function is
/// passed on to the ModelInitializer creating random positions according to
/// some distribution.
TumorCell* CreateTumorCell(const Double3& position);

/// Wrapper to multiple call to CreateTumorCell.
void PlaceTumorCells(std::vector<Double3>& positions);

/// Compute the volume of a cylinder with given diameter and length
double CylinderVolume(double diameter, double length);

/// @brief  This function places a vessel in the simulation
/// @param start Beginning of the vessel
/// @param end end of the vessel
/// @param compartment_length Length of the individual compartments (agents)
double PlaceStraightVessel(
    Double3 start, Double3 end, double compartment_length, double diameter,
    AgentPointer<neuroscience::NeuronOrNeurite>* parent = nullptr,
    AgentPointer<neuroscience::NeuronOrNeurite>* terminal_vessel = nullptr);

/// @brief This function places a random vessel in the simulation. Uses
/// parameters from the simulation parameters.
/// @param start Beginning of the vessel
/// @param end end of the vessel
/// @param param Simulation parameters
double PlaceRandomVessel(Double3 start, Double3 end, double diameter,
                         unsigned int random_seed);

/// @brief 3D Gaussian function
double Gaussian(double x, double y, double z);

/// @brief  This function sets up the experiments.
/// @param experiment experiment to be set up
/// @param fn Initial concentration of nutrients
/// @param fv Initial concentration of VEGF
/// @param fd Initial concentration of DOX
/// @param ft Initial concentration of TRA
/// @param bct_n Boundary condition type for nutrients
/// @param bct_v Boundary condition type for VEGF
/// @param bct_d Boundary condition type for DOX
/// @param bct_t Boundary condition type for TRA
/// @param initialize_random_cells Whether to initialize cells randomly
/// @param initialize_tumor_spheroid Whether to initialize tumor spheroid
/// @param initialize_vasculature Whether to initialize vasculature
/// @param sparam
void SetUpExperiment(const Experiment experiment,
                     std::function<double(double, double, double)>& fn,
                     std::function<double(double, double, double)>& fv,
                     std::function<double(double, double, double)>& fd,
                     std::function<double(double, double, double)>& ft,
                     BoundaryConditionType& bct_n, BoundaryConditionType& bct_v,
                     BoundaryConditionType& bct_d, BoundaryConditionType& bct_t,
                     bool& initialize_random_cells,
                     bool& initialize_tumor_spheroid,
                     bool& initialize_vasculature, const Param* param,
                     const SimParam* sparam);

/// @brief  This function initializes the vessel structure in the simulation
/// @param experiment experiment to be set up
/// @param sparam Simulation parameters
void InitializeVessels(const Experiment experiment, const Param* param,
                       const SimParam* sparam);

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

  constexpr Experiment experiment = Experiment::kFullScaleModel;

  // ---------------------------------------------------------------------------
  // 2. Define parameters and initialize simulation
  // ---------------------------------------------------------------------------
  auto set_param = [](Param* param) {
    param->calculate_gradients = true;
    param->visualization_interval = static_cast<uint32_t>(std::floor(
        static_cast<double>(param->Get<SimParam>()->visualization_interval) /
        param->simulation_time_step));
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
  // 3. Get setup for experiment
  // ---------------------------------------------------------------------------

  // Initial concentrations of the different substances
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

  // Boundary conditions for the substances
  auto bc_nutrients = std::make_unique<ConstantBoundaryCondition>(
      sparam->boundary_condition_nutrients);
  auto bc_vegf = std::make_unique<ConstantBoundaryCondition>(
      sparam->boundary_condition_vegf);
  auto bc_dox = std::make_unique<ConstantBoundaryCondition>(
      sparam->boundary_condition_dox);
  auto bc_tra = std::make_unique<ConstantBoundaryCondition>(
      sparam->boundary_condition_tra);

  // Boudary condition types for the substances
  auto bc_type_nutrients = BoundaryConditionType::kNeumann;
  auto bc_type_vegf = BoundaryConditionType::kNeumann;
  auto bc_type_dox = BoundaryConditionType::kNeumann;
  auto bc_type_tra = BoundaryConditionType::kNeumann;

  // Initialize cells randomly at the beginning of the simulation
  bool initialize_random_cells = false;
  // Initialize tumor spheroid at the beginning of the simulation
  bool initialize_tumor_spheroid = false;
  // Initialize vasculature at the beginning of the simulation
  bool initialize_vasculature = false;

  SetUpExperiment(
      experiment, initial_nutrient_concentration, initial_vegf_concentration,
      initial_dox_concentration, initial_tra_concentration, bc_type_nutrients,
      bc_type_vegf, bc_type_dox, bc_type_tra, initialize_random_cells,
      initialize_tumor_spheroid, initialize_vasculature, param, sparam);

  // ---------------------------------------------------------------------------
  // 4. Define continuum models for nutrients and VEGF
  // ---------------------------------------------------------------------------

  // Define nutrients with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kNutrients, "Nutrients", sparam->diffusion_nutrients,
      sparam->decay_rate_nutrients, sparam->diffusion_resolution_nutrients);
  ModelInitializer::InitializeSubstance(Substances::kNutrients,
                                        initial_nutrient_concentration);
  ModelInitializer::AddBoundaryConditions(
      Substances::kNutrients, bc_type_nutrients, std::move(bc_nutrients));

  // Define VEGF with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kVEGF, "VEGF", sparam->diffusion_vegf,
      sparam->decay_rate_vegf, sparam->diffusion_resolution_vegf);
  ModelInitializer::InitializeSubstance(Substances::kVEGF,
                                        initial_vegf_concentration);
  ModelInitializer::AddBoundaryConditions(Substances::kVEGF, bc_type_vegf,
                                          std::move(bc_vegf));

  // Define TRA with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kTRA, "TRA", sparam->diffusion_tra, sparam->decay_rate_tra,
      sparam->diffusion_resolution_tra);
  ModelInitializer::InitializeSubstance(Substances::kTRA,
                                        initial_tra_concentration);
  ModelInitializer::AddBoundaryConditions(Substances::kTRA, bc_type_tra,
                                          std::move(bc_tra));

  // Define DOX with constant initial conditions
  ModelInitializer::DefineSubstance(
      Substances::kDOX, "DOX", sparam->diffusion_dox, sparam->decay_rate_dox,
      sparam->diffusion_resolution_dox);
  ModelInitializer::InitializeSubstance(Substances::kDOX,
                                        initial_dox_concentration);
  ModelInitializer::AddBoundaryConditions(Substances::kDOX, bc_type_dox,
                                          std::move(bc_dox));

  // Define upper and lower threshold for nutrients
  rm->ForEachDiffusionGrid([](DiffusionGrid* grid) {
    grid->SetUpperThreshold(1.0);
    grid->SetLowerThreshold(0.0);
  });

  // Deactivate gradient calculation for nutrients, TRA, and DOX
  rm->GetDiffusionGrid(Substances::kNutrients)->TurnOffGradientCalculation();
  rm->GetDiffusionGrid(Substances::kTRA)->TurnOffGradientCalculation();
  rm->GetDiffusionGrid(Substances::kDOX)->TurnOffGradientCalculation();

  // ---------------------------------------------------------------------------
  // 5. Define initial configurations of agents
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
      ModelInitializer::CreateAgentsInSphereRndm(
          {0, 0, 0}, sparam->GetSpheroidRadius(), sparam->num_cells,
          CreateTumorCell);
    }

    if (initialize_vasculature) {
      InitializeVessels(experiment, param, sparam);
    }
  }

  // ---------------------------------------------------------------------------
  // 6. Track simulation results over time with timeseries objects
  // ---------------------------------------------------------------------------

  // Collect the number of Cells in different states over time
  DefineAndRegisterCollectors();

  // ---------------------------------------------------------------------------
  // 7. Use force module typically used by UT Austin
  // ---------------------------------------------------------------------------

  // Use custom force module implemented in MechanicalInteractionForce
  // Note that the force module currently does not support any forces
  // between vessels and cells.
  auto* custom_force = new MechanicalInteractionForce(
      sparam->adhesion_scale_parameter, sparam->repulsive_scale_parameter);
  auto* op = scheduler->GetOps("mechanical forces")[0];
  op->frequency_ = sparam->force_calculation_frequency;
  auto* force_implementation = op->GetImplementation<MechanicalForcesOp>();
  force_implementation->SetInteractionForce(custom_force);

  // ---------------------------------------------------------------------------
  // 8. Specific fix for force and environment combination
  // ---------------------------------------------------------------------------

  // Set box length manually because our interaction range is larger than the
  // cell's diameter. In the current setup we restrict vessel growth once we
  // come close to a tumor cell. We set the UniformGrid to larger box sizes
  // such that we can resolve more long distance relationships than with the
  // UniformGrid.
  double distance_for_growth_stop = 60;
  double box_length =
      std::ceil(2 * sparam->action_radius_factor * sparam->cell_radius);
  box_length = std::max(box_length, distance_for_growth_stop);
  env->SetBoxLength(static_cast<int32_t>(box_length));

  // ---------------------------------------------------------------------------
  // 9. Load balance (Linux only)
  // ---------------------------------------------------------------------------
#ifdef __linux__
  scheduler->GetOps("load balancing")[0]->frequency_ = 20;
#endif  // __linux__

  // ---------------------------------------------------------------------------
  // 10. Track continuum models
  // ---------------------------------------------------------------------------

  if (sparam->verify_continuum_values) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "VerifyContinuum", OpComputeTarget::kCpu, new VerifyContinuum());
    auto* verify_continuum = NewOperation("VerifyContinuum");
    scheduler->ScheduleOp(verify_continuum, OpType::kPostSchedule);
  }

  // ---------------------------------------------------------------------------
  // 11. Track continuum models
  // ---------------------------------------------------------------------------

  if (sparam->tip_cell_finder_update_frequency <
      std::numeric_limits<int>::max()) {
    OperationRegistry::GetInstance()->AddOperationImpl(
        "update tip-cell finder", OpComputeTarget::kCpu,
        new UpdateTipCellFinder());
    auto* update_tip_cell_finder = NewOperation("update tip-cell finder");
    update_tip_cell_finder->frequency_ =
        sparam->tip_cell_finder_update_frequency;
    scheduler->ScheduleOp(update_tip_cell_finder, OpType::kPostSchedule);
  }

  // ---------------------------------------------------------------------------
  // 12. Run simulation and visualize results
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

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

TumorCell* CreateTumorCell(const Double3& position) {
  // Connect to active simulation, get parameters and a random generator
  auto* sim = Simulation::GetActive();
  auto* param = sim->GetParam();
  auto* sparam = param->Get<SimParam>();
  // Create cells at a random position (randomness through model initializer)
  // that can divide and are quiescent
  int cell_state = CellState::kQuiescent;
  auto* tumor_cell = new TumorCell(position, cell_state);
  // Set radius, nuclear radius, and action radius
  tumor_cell->SetActionRadiusFactor(sparam->action_radius_factor);
  tumor_cell->SetRadii(sparam->cell_radius, sparam->cell_nuclear_radius,
                       sparam->action_radius_factor * sparam->cell_radius);
  // Cells gain half their volume during the growth phase.
  double growth_rate = 2.0 / 3.0 * Math::kPi * pow(sparam->cell_radius, 3) /
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
  for (const auto& pos : positions) {
    auto* tumor_cell = CreateTumorCell(pos);
    rm->AddAgent(tumor_cell);
  }
}

double CylinderVolume(double diameter, double length) {
  return Math::kPi * diameter * diameter * length / 4.0;
}

double PlaceStraightVessel(
    Double3 start, Double3 end, double compartment_length, double diameter,
    AgentPointer<neuroscience::NeuronOrNeurite>* parent,
    AgentPointer<neuroscience::NeuronOrNeurite>* terminal_vessel) {
  auto* rm = Simulation::GetActive()->GetResourceManager();
  auto* param = Simulation::GetActive()->GetParam();
  auto* sparam = param->Get<SimParam>();

  // Compute parameters for straight line between start and end.
  Double3 direction = end - start;
  double distance = direction.Norm();
  direction.Normalize();
  auto n_compartments =
      static_cast<int>(std::floor(distance / compartment_length));
  compartment_length = distance / n_compartments;
  // double delta = distance - n_compartments * compartment_length;
  double delta = 0.0;
  std::cout << "delta: " << delta << std::endl;

  double vessel_volume = CylinderVolume(diameter, distance);

  // // Warn if chosen parameters are not selected ideally
  // if (abs(n_compartments * compartment_length - distance) > 1e-2) {
  //   Log::Warning("PlaceStraightVessel",
  //                "Vessel will be shorter than expected.");
  // }

  // The setup requires us to define a NeuronSoma, which is kind of a left over
  // from the neuroscience module.

  // Define a first neurite
  Vessel v;  // Used for prototype argument (virtual+template not supported c++)
  Vessel* vessel_compartment_1{nullptr};
  double connection_factor = 0.5;
  bool start_new_vessel = true;
  if (parent != nullptr && *parent != nullptr) {
    std::cout << "Connecting to parent" << std::endl;
    vessel_compartment_1 = new Vessel();
    vessel_compartment_1->SetMother(*parent);
    rm->AddAgent(vessel_compartment_1);
    auto* parent_vessel = dynamic_cast<Vessel*>((*parent).Get());
    if (parent_vessel->GetDaughterLeft() == nullptr) {
      parent_vessel->SetDaughterLeft(
          vessel_compartment_1->GetAgentPtr<neuroscience::NeuriteElement>());
    } else if (parent_vessel->GetDaughterRight() == nullptr) {
      parent_vessel->SetDaughterRight(
          vessel_compartment_1->GetAgentPtr<neuroscience::NeuriteElement>());
    } else {
      Log::Warning("PlaceStraightVessel",
                   "Parent vessel already has two daughters.");
      start_new_vessel = true;
    }
    vessel_compartment_1->SetRestingLength(compartment_length);
    vessel_compartment_1->SetSpringAxis(direction);
    connection_factor = 0.5;
    start_new_vessel = false;
  }
  if (start_new_vessel) {
    std::cout << "Starting new vessel" << std::endl;
    const Double3 tmp = start;
    auto* soma = new neuroscience::NeuronSoma(tmp);
    rm->AddAgent(soma);
    vessel_compartment_1 =
        dynamic_cast<Vessel*>(soma->ExtendNewNeurite(direction, &v));
  }
  // vessel_compartment_1->SetPosition(start + direction * compartment_length *
  //                                               connection_factor);
  vessel_compartment_1->SetPosition(start);
  vessel_compartment_1->SetMassLocation(start + direction * compartment_length);
  vessel_compartment_1->SetActualLength(2 * (1 - connection_factor) *
                                        compartment_length);
  vessel_compartment_1->SetDiameter(diameter);
  vessel_compartment_1->ProhibitGrowth();

  Vessel* vessel_compartment_2{nullptr};
  for (int i = 1; i < n_compartments; i++) {
    double extension = 0;
    if (i == n_compartments - 1) {
      extension = delta / compartment_length;
    }
    // Compute location of next vessel element
    Double3 agent_position =
        start + direction * compartment_length * (static_cast<double>(i) + 0.5);
    Double3 agent_end_position =
        start + direction * compartment_length *
                    (static_cast<double>(i) + 1.0 + extension);
    // Create new vessel
    vessel_compartment_2 = new Vessel();
    // Set position an length
    vessel_compartment_2->SetPosition(agent_position);
    vessel_compartment_2->SetMassLocation(agent_end_position);
    vessel_compartment_2->SetActualLength(compartment_length *
                                          (1.0 + extension));
    vessel_compartment_2->SetRestingLength(compartment_length *
                                           (1.0 + extension));
    vessel_compartment_2->SetSpringAxis(direction);
    vessel_compartment_2->SetDiameter(diameter);
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
  if (terminal_vessel != nullptr) {
    *terminal_vessel =
        vessel_compartment_1->GetAgentPtr<neuroscience::NeuronOrNeurite>();
  }
  return vessel_volume;
}

double PlaceRandomVessel(Double3 start, Double3 end, double diameter,
                         unsigned int random_seed) {
  // Get parameters
  auto* rm = Simulation::GetActive()->GetResourceManager();
  auto* param = Simulation::GetActive()->GetParam();
  auto* sparam = param->Get<SimParam>();
  auto* nparam = param->Get<neuroscience::Param>();

  // Track vessel volume
  double vessel_volume = 0.0;

  // Compute parameters for straight line between start and end.
  Double3 global_direction = end - start;
  double distance = global_direction.Norm();
  global_direction.Normalize();

  // Get an orthogonal system to the direction vector
  Double3 ortho1;
  Double3 ortho2;
  GetOrthogonalSystem(global_direction, ortho1, ortho2);

  // Get two random fields realizations for the vessel
  RandomField rf(sparam->random_vessel_num_modes, distance,
                 2 * nparam->neurite_min_length, sparam->random_vessel_exponent,
                 sparam->random_vessel_max_deviation * distance,
                 sparam->random_vessel_frequency_mean,
                 sparam->random_vessel_frequency_std, random_seed);
  std::vector<double> random_field_1;
  std::vector<double> random_field_2;
  rf.GetRealization(random_field_1);
  rf.GetRealization(random_field_2);

  // Compute the discretization length along the straight line
  const double discretization_length = distance / rf.GetNumPoints();

  // The setup requires us to define a NeuronSoma, which is kind of a left over
  // from the neuroscience module.
  const Double3 tmp = start;
  auto* soma = new neuroscience::NeuronSoma(tmp);
  rm->AddAgent(soma);

  // Define positions and lengths of the agent
  Double3 agent_position_start = start;
  Double3 offset = ortho1 * random_field_1[1] + ortho2 * random_field_2[1];
  Double3 agent_position_end =
      start + global_direction * discretization_length + offset;
  Double3 agent_direction = agent_position_end - agent_position_start;
  double compartment_length = agent_direction.Norm();

  // Track vessel volume
  vessel_volume += CylinderVolume(compartment_length, diameter);

  // Define a first neurite
  Vessel v;  // Used for prototype argument (virtual+template not supported c++)
  auto* vessel_compartment_1 =
      dynamic_cast<Vessel*>(soma->ExtendNewNeurite(agent_direction, &v));
  vessel_compartment_1->SetPosition(start +
                                    agent_direction * compartment_length * 0.5);
  vessel_compartment_1->SetMassLocation(start +
                                        agent_direction * compartment_length);
  vessel_compartment_1->SetActualLength(compartment_length);
  vessel_compartment_1->SetDiameter(diameter);
  vessel_compartment_1->ProhibitGrowth();

  Vessel* vessel_compartment_2{nullptr};
  for (int i = 1; i < rf.GetNumPoints() - 2; i++) {
    // Compute location of next vessel element
    agent_position_start = agent_position_end;
    offset = ortho1 * random_field_1[i + 1] + ortho2 * random_field_2[i + 1];
    agent_position_end =
        start + global_direction * (i + 1) * discretization_length + offset;
    agent_direction = agent_position_end - agent_position_start;
    compartment_length = agent_direction.Norm();

    // Track vessel volume
    vessel_volume += CylinderVolume(compartment_length, diameter);

    Double3 agent_position =
        agent_position_start + agent_direction * compartment_length * 0.5;
    Double3 agent_mass_position = agent_position_end;
    // Create new vessel
    vessel_compartment_2 = new Vessel();
    // Set position an length
    vessel_compartment_2->SetPosition(agent_position);
    vessel_compartment_2->SetMassLocation(agent_mass_position);
    vessel_compartment_2->SetActualLength(compartment_length);
    vessel_compartment_2->SetRestingLength(compartment_length);
    vessel_compartment_2->SetSpringAxis(agent_direction);
    vessel_compartment_2->SetDiameter(diameter);
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
  return vessel_volume;
};

double Gaussian(double x, double y, double z) {
  double mu_x = 0.0;
  double mu_y = 0.0;
  double mu_z = 0.0;
  double sigma = 130.0;
  double r = std::sqrt(std::pow(x - mu_x, 2) + std::pow(y - mu_y, 2) +
                       std::pow(z - mu_z, 2));

  return std::exp(-r * r / (2 * sigma * sigma));
}

void SetUpExperiment(const Experiment experiment,
                     std::function<double(double, double, double)>& fn,
                     std::function<double(double, double, double)>& fv,
                     std::function<double(double, double, double)>& fd,
                     std::function<double(double, double, double)>& ft,
                     BoundaryConditionType& bct_n, BoundaryConditionType& bct_v,
                     BoundaryConditionType& bct_d, BoundaryConditionType& bct_t,
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
      bct_n = BoundaryConditionType::kDirichlet;
      bct_v = BoundaryConditionType::kNeumann;
      bct_d = BoundaryConditionType::kOpenBoundaries;
      bct_t = BoundaryConditionType::kOpenBoundaries;
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
      bct_n = BoundaryConditionType::kDirichlet;
      bct_v = BoundaryConditionType::kNeumann;
      bct_d = BoundaryConditionType::kOpenBoundaries;
      bct_t = BoundaryConditionType::kOpenBoundaries;
      initialize_tumor_spheroid = true;
      initialize_vasculature = false;
      break;
    case Experiment::kSpheroidTreatment:
      bct_n = BoundaryConditionType::kClosedBoundaries;
      bct_v = BoundaryConditionType::kClosedBoundaries;
      bct_v = BoundaryConditionType::kClosedBoundaries;
      bct_d = BoundaryConditionType::kClosedBoundaries;
      initialize_tumor_spheroid = true;
      initialize_vasculature = false;
      break;
    case Experiment::kVesselsToCenter:
      fv = Gaussian;
      bct_n = BoundaryConditionType::kOpenBoundaries;
      bct_v = BoundaryConditionType::kDirichlet;
      bct_v = BoundaryConditionType::kOpenBoundaries;
      bct_d = BoundaryConditionType::kOpenBoundaries;
      initialize_tumor_spheroid = false;
      initialize_vasculature = true;
      break;
    case Experiment::kVesselsCoupling:
      fv = [slope, offset](double x, double, double) {
        return slope * x + offset;
      };
      bct_n = BoundaryConditionType::kNeumann;
      bct_v = BoundaryConditionType::kOpenBoundaries;
      bct_v = BoundaryConditionType::kOpenBoundaries;
      bct_d = BoundaryConditionType::kOpenBoundaries;
      initialize_tumor_spheroid = false;
      initialize_vasculature = true;
      break;
    case Experiment::kSimplifiedGrowth:
      Log::Fatal("SetUpExperiment", "Not implemented yet");
      break;
    case Experiment::kFullScaleModel:
      bct_n = BoundaryConditionType::kNeumann;
      bct_v = BoundaryConditionType::kNeumann;
      bct_v = BoundaryConditionType::kNeumann;
      bct_d = BoundaryConditionType::kNeumann;
      initialize_tumor_spheroid = true;
      initialize_vasculature = true;
      break;
    default:
      Log::Fatal("SetUpExperiment", "Unknown experiment");
  };

  // Check if all functions are set
  if (!fn || !fv || !fd || !ft) {
    Log::Fatal("SetUpExperiment", "Not all functions are set");
  }
};

void InitializeVessels(const Experiment experiment, const Param* param,
                       const SimParam* sparam) {
  if (experiment == Experiment::kVesselsToCenter) {
    PlaceStraightVessel({-200, 0, -400}, {-200, 0, 400},
                        sparam->default_vessel_length, 15);
    PlaceStraightVessel({200, 0, -400}, {200, 0, 400},
                        sparam->default_vessel_length, 15);
    PlaceStraightVessel({0, -400, 200}, {0, 400, 200},
                        sparam->default_vessel_length, 15);
    PlaceStraightVessel({0, -400, -200}, {0, 400, -200},
                        sparam->default_vessel_length, 15);
  } else if (experiment == Experiment::kVesselsCoupling) {
    PlaceStraightVessel({0, 0, -400}, {0, 0, 400},
                        sparam->default_vessel_length, 15);
  } else if (experiment == Experiment::kFullScaleModel) {
    //   auto* rnd = Simulation::GetActive()->GetRandom();

    //   // Parameters extracted from rat brain (Secomb)
    //   constexpr double vessels_per_volume = 104.0;
    //   constexpr double volume = 550 * 550 * 230;  // in um^3
    //   constexpr double gev_location = 7.53401234807571;
    //   constexpr double gev_scale = 2.3153691046974907;
    //   constexpr double gev_xi = -0.3292535751714517;
    //   constexpr double gev_min = 5.0;   // for numerical integration
    //   constexpr double gev_max = 50.0;  // for numerical integration
    //   constexpr double wald_location = 10.977589465183868;
    //   constexpr double wald_scale = 63.81303941790211;
    //   constexpr double wald_min = 100;     // for numerical integration
    //   constexpr double wald_max = 1000.0;  // for numerical integration
    //   constexpr double random_vessel_threshold = 200;

    //   // Compute the fraction of the vessels that we want to neglect (because
    //   to
    //   // small). Typical numbers: wald_min = 50 -> 45%, wald_min = 100 ->
    //   78%.
    //   // Use Gauss-Kronrod 21-point integration rule.
    //   ROOT::Math::IntegratorOneDim integrator;
    //   auto g = [](double x) { return wald_pdf(x, wald_location, wald_scale);
    //   }; double wald_integral_to_min = integrator.Integral(g, 0.0, wald_min);

    //   // Parameters
    //   const auto min = param->min_bound;
    //   const auto max = param->max_bound;
    //   const double simulation_volume = std::pow(max - min, 3);

    //   // Compute number of vessels
    //   const int num_vessels = static_cast<int>(
    //       std::ceil(vessels_per_volume * (1 - wald_integral_to_min) *
    //                 simulation_volume / volume));

    //   // Print computed info
    //   std::cout << "Number of vessels: " << num_vessels << std::endl;
    //   std::cout << "Simulation volume: " << simulation_volume << std::endl;
    //   std::cout << "Volume: " << volume << std::endl;
    //   std::cout << "Min: " << min << std::endl;
    //   std::cout << "Max: " << max << std::endl;
    //   std::cout << "Wald integral to min: " << wald_integral_to_min <<
    //   std::endl;

    //   // Generate random vessel diameters, lengths, and starting points
    //   std::vector<double> vessel_diameters(num_vessels);
    //   std::vector<double> vessel_lengths(num_vessels);
    //   std::vector<Double3> vessel_start_points(num_vessels);
    //   std::vector<Double3> vessel_end_points(num_vessels);

    //   // Track vessel volume
    //   double vessel_volume = 0.0;

    //   // Generate random vessel diameters via the GENEXTREME distribution
    //   auto gev_distribution = [](const double* x, const double* params) {
    //     return gev_pdf(x[0], params[0], params[1], params[2]);
    //   };
    //   auto rng_gev = rnd->GetUserDefinedDistRng1D(
    //       gev_distribution, {gev_location, gev_scale, gev_xi}, gev_min,
    //       gev_max);
    //   for (int i = 0; i < num_vessels; i++) {
    //     vessel_diameters[i] = rng_gev.Sample();
    //   }

    //   // Generate random vessel lengths via the WALD distribution
    //   auto wald_distribution = [](const double* x, const double* params) {
    //     return wald_pdf(x[0], params[0], params[1]);
    //   };
    //   auto rng_wald = rnd->GetUserDefinedDistRng1D(
    //       wald_distribution, {wald_location, wald_scale}, wald_min,
    //       wald_max);
    //   for (int i = 0; i < num_vessels; i++) {
    //     vessel_lengths[i] = rng_wald.Sample();
    //   }

    //   // Generate random vessel starting points
    //   for (int i = 0; i < num_vessels; i++) {
    //     vessel_start_points[i] = {rnd->Uniform(min, max), rnd->Uniform(min,
    //     max),
    //                               rnd->Uniform(min, max)};
    //     // Generate random vessel end points
    //     vessel_end_points[i] =
    //         rnd->Sphere(vessel_lengths[i]) + vessel_start_points[i];
    //     while (vessel_end_points[i][0] < min || vessel_end_points[i][0] > max
    //     ||
    //            vessel_end_points[i][1] < min || vessel_end_points[i][1] > max
    //            || vessel_end_points[i][2] < min || vessel_end_points[i][2] >
    //            max) {
    //       vessel_end_points[i] =
    //           rnd->Sphere(vessel_lengths[i]) + vessel_start_points[i];
    //     }
    //   }

    //   // Plot histogram of the diameters and lengths
    //   std::string diameters_file = "diameters";
    //   std::string lengths_file = "lengths";
    //   PlotAndSaveHistogram(vessel_diameters, diameters_file);
    //   PlotAndSaveHistogram(vessel_lengths, lengths_file);

    //   // Define random seed
    //   unsigned int seed = 0;

    //   for (int i = 0; i < num_vessels; i++) {
    //     const auto& start = vessel_start_points[i];
    //     const auto& end = vessel_end_points[i];
    //     const auto& length = vessel_lengths[i];
    //     const auto& diameter = vessel_diameters[i];
    //     // const auto diameter = 15;
    //     if (length < random_vessel_threshold) {
    //       vessel_volume += PlaceStraightVessel(
    //           start, end, sparam->default_vessel_length, diameter);
    //     } else {
    //       vessel_volume += PlaceRandomVessel(start, end, diameter, seed++);
    //     }
    //   }
    //   // Print vessel volumen and vessel volume fraction
    //   std::cout << "Vessel volume: " << vessel_volume << std::endl;
    //   std::cout << "Vessel volume fraction: " << vessel_volume /
    //   simulation_volume
    //             << std::endl;

    constexpr bool kWithConnectivity = true;
    DataParserVTP parser;
    parser.ParseData("data/network.vtp");
    // parser.ParseData("data/dummy.vtp");
    parser.PostProcessData();

    if (kWithConnectivity) {
      // ---------------------------------------------------------------------
      // With connectivity
      // ---------------------------------------------------------------------
      const auto& num_lines = parser.GetNumLines();
      const auto& points = parser.GetPoints();
      const auto& connectivity = parser.GetConnectivity();
      const auto& radius = parser.GetRadii();
      std::vector<AgentPointer<neuroscience::NeuronOrNeurite>> terminal_ends(
          points.size());
      for (size_t i = 0; i < num_lines; i++) {
        // Define start and end points
        size_t start_index = connectivity[2 * i];
        size_t end_index = connectivity[2 * i + 1];
        if (terminal_ends[start_index] == nullptr &&
            terminal_ends[end_index] != nullptr) {
          std::swap(start_index, end_index);
        }
        const auto& start = points[start_index];
        const auto& end = points[end_index];
        const auto& diameter = 2 * radius[i];
        // PlaceStraightVessel(start, end, sparam->default_vessel_length,
        // diameter);
        // PlaceStraightVessel(start, end, sparam->default_vessel_length,
        // diameter,
        //                     &terminal_ends[start_index],
        //                     &terminal_ends[end_index]);
        PlaceStraightVessel(start, end, 2.0, diameter,
                            &terminal_ends[start_index],
                            &terminal_ends[end_index]);
      }
    } else {
      // ---------------------------------------------------------------------
      // Without connectivity
      // ---------------------------------------------------------------------
      for (auto& segment : parser.data) {
        const auto& start = segment.start_position;
        const auto& end = segment.end_position;
        const auto& diameter = 2 * segment.radius;
        // PlaceStraightVessel(start, end, sparam->default_vessel_length,
        //                     diameter);
        PlaceStraightVessel(start, end, 2.0, diameter);
      }
    }
  } else {
    Log::Fatal("InitializeVessels",
               "No vessel structure defined for this "
               "experiment");
  }
}

}  // namespace bdm
