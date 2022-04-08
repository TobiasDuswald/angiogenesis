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

// This file collect the parameters needed for the simulation. The parameters
// are mostly taken form Rocha et al 2018 and Lima et al 2021. The majority of
// the parameter belong to formation of the tumor which is not considered in the
// early stages of the project.

#ifndef SIM_PARAM_H_
#define SIM_PARAM_H_

#include "biodynamo.h"

namespace bdm {

// Available substances in Simulation
enum Substances { kNutrients, kVEGF };

// This class defines parameters that are specific to this simulation. The unit
// h refers to hours.
struct SimParam : public ParamGroup {
  BDM_PARAM_GROUP_HEADER(SimParam, 1);

  // Total simulation time (unit [min]). This unit carries over to
  // bdm::Param.simulation_time_step. E.g. a timestep of 0.01 [min] = 0.6 sec.
  // The parameters are chosen such that no cell can move more than 0.2 \mu m
  // per simulation step for typical forces. For the simulation time below,
  // specify the days in  <days>*<24 hours>*<60 minutes>.
  double total_sim_time{0.01 * 24 * 60};

  // If running on linux, this parameter determines how often we export a
  // paraview visualization, e.g. the interval between exports in unit [min].
  double visualization_interval{1.0};

  // // Number of (Tumor)Cells that we create in the domain at the beginning of
  // // the simulation (uint [1])
  // u_int64_t no_cells{1000};

  // Length of vessels at initialization
  double default_vessel_length{10};

  // Lower bound for the domain (applies to x,y,z; unit [\mu m])
  double lower_bound{-200.0};

  // Upper bound for the domain (applies to x,y,z; unit [\mu m])
  double upper_bound{200.0};

  // Cell radius (unit [\mu m])
  double cell_radius{9.953};

  // Cell radius standard deviation for random gaussian initialization
  // (unit [\mu m])
  double cell_radius_sigma{0.4};

  // Radius of nucleus (unit [\mu m])
  double cell_nuclear_radius{5.296};

  // Factor to compute action radius from actual cell radius. R_A = \alpha R.
  // (unit [1])
  double action_radius_factor{1.214};

  // Duration of the cell cycle (numeric parameter \tau_{P}, unit [min])
  double duration_cell_cycle{18.0 * 60.0};

  // Duration of the growth phase (numeric parameter \tau_{G1}, unit [h])
  double duration_growth_phase{9.0 * 60.0};

  // Duration of apoptosis
  double duration_apoptosis{8.6 * 60.0};

  // Apoptosis rate (numeric parameter \bar{\alpha}_{D}, unit [min^{-1}])
  double apoptosis_rate{0.000408 / 60.0};

  // Undefined rate governing the transition from quiescent to proliferative
  // (numeric parameter \bar{\alpha}_{P}, unit [min^{-1}])
  double qp_transition_rate{0.0493 / 60.0};

  // Hypoxic threshold (numeric parameter \sigma_{H}, unit [1])
  double hypoxic_threshold{0.0538};

  // Mathematical parameter k (unit [1])
  double k{50.0};

  // Gamma factor relating cell death and lack of nutrients (unit [min^{-1}]))
  double gamma{0.0245 / 60.0};

  // Viscosity of the surrounding (numerical parameter \nu - not defined in
  // paper, taken/assumed from Git)
  double viscosity{2.0};

  // Maximum speed that cells can move with (not defined in
  // paper, seen in Git, unit [\mu m / min]). With this maximum velocity and
  // the BDM default timestep 0.01, cells move at most 0.1 \mu m per timestep.
  double max_speed{10.0};

  // Numeric parameter c_{cca} for force, unit [\mu m / min].
  double adhesion_scale_parameter{0.0489};

  // Numeric parameter c_{ccr}for force, unit [\mu m / min].
  double repulsive_scale_parameter{10.0};

  // Resolution of the diffusion grid. Gets forwarded to DiffusionGrid
  // constructor
  int diffusion_resolution{50};

  // Diffusion coefficient for glucose (unit [\mu m / min])
  double diffusion_nutrients{50.0 / 60.0};

  // Uptake rate of glucose by cells (unit [min^{-1}])
  double uptake_rate_glucose{0.0483 / 60.0};

  // Initial value of the glucose concentarion, uniform over grid (unit [1])
  double initial_nutrient_concentration{0.5};

  // The decay constant of glucose. It basically describes an exponential
  // decay for each point in the diffusion grid.
  double decay_rate_nutrients{0.00001};

  // Secretion rate tumor cells, i.e. how much VEGF is released by a tumor cell
  // per minute.
  double secretion_rate_vegf{0.03 / 60.0};

  // Diffusion constant for VEGF
  double diffusion_vegf{40.0 / 60.0};

  // Decay constant for VEGF
  double decay_rate_vegf{0.0};

  // Parameter to decide if dead cells decrease in size and are removed or if we
  // keep them in the simulation.
  bool keep_dead_cells{false};

  //////////////////////////////////////////////////////////////////////////////
  // Vessel parameters
  //////////////////////////////////////////////////////////////////////////////

  // VEGF threshold for sprouting
  double vegf_threshold_sprouting{1e-3};

  // VEGF gradient threshold for apical growth
  double vegf_grad_threshold_apical_growth{1e-5};

  // Minimum distance to bifurcatoin or terminal end of vessel to allow
  // sprouting
  double min_dist_to_bifurcation{60.0};

  // Sprouting probability
  double sprouting_probability{0.001};

  // Weight for random direction of the apical growth
  double apical_growth_random_weight{0.2};

  // Weight for old direction of the apical growth
  double apical_growth_old_weight{0.5};

  // Weight for gradient direction of the apical growth
  double apical_growth_gradient_weight{0.3};

  // Apical growth speed
  double apical_growth_speed{1.0};

  // Nutrient supply by vessel (unit Nutrients / (Area * min)])
  double nutrient_supply_rate_vessel{0.0001};
};

}  // namespace bdm

#endif  // SIM_PARAM_H_
