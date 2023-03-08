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
enum Substances { kNutrients, kVEGF, kTRA, kDOX };

// Different setups for the computational experiments.
enum class Experiment {
  kAvascularTumorSpheroid,
  kPorousTumorSpheroid,
  kSpheroidTreatment,
  kVesselsToCenter,
  kVesselsCoupling,
  kSimplifiedGrowth,
  kFullScaleModel
};

// This class defines parameters that are specific to this simulation. The unit
// h refers to hours.
struct SimParam : public ParamGroup {
  BDM_PARAM_GROUP_HEADER(SimParam, 1);

  // -----------------------------------------------------------------------
  // Simulation parameters
  // -----------------------------------------------------------------------

  // Total simulation time (unit [min]). This unit carries over to
  // bdm::Param.simulation_time_step. E.g. a timestep of 0.01 [min] = 0.6 sec.
  // The parameters are chosen such that no cell can move more than 0.2 \mu m
  // per simulation step for typical forces. For the simulation time below,
  // specify the days in  <days>*<24 hours>*<60 minutes>.
  double total_sim_time{0.01 * 24 * 60};

  // This parameter determines how often we export a paraview visualization,
  // e.g. the interval between exports in unit [min]. Used to set
  // bdm::Param::visualization_interval.
  double visualization_interval{1.0};

  // Number of (Tumor)Cells that we create in the domain at the beginning of
  // the simulation (uint [1])
  u_int64_t num_cells{1000};

  // Volume filling for the tumor spheroid (unit [1]). Used to determine the
  // initial radius of the tumor spheroid.
  double filled_volume{0.7};

  // Parameter to decide if dead cells decrease in size and are removed or if we
  // keep them in the simulation.
  bool keep_dead_cells{true};

  // Verify that the continuum values are all between 0 and 1
  bool verify_continuum_values{true};

  // Update frequency for TipCellFinder
  size_t tip_cell_finder_update_frequency{1};

  // -----------------------------------------------------------------------
  // TumorCell parameters
  // -----------------------------------------------------------------------

  // Cell radius (unit [\mu m])
  double cell_radius{9.953};

  // Radius of nucleus (unit [\mu m])
  double cell_nuclear_radius{5.296};

  // Factor to compute action radius from actual cell radius. R_A = \alpha R.
  // (unit [1])
  double action_radius_factor{1.214};

  // -----------------------------------------------------------------------
  // Cell cycle parameters
  // -----------------------------------------------------------------------

  // Duration of the cell cycle (numeric parameter \tau_{P}, unit [min])
  double duration_cell_cycle{18.0 * 60.0};

  // Duration of the growth phase (numeric parameter \tau_{G1}, unit [h])
  double duration_growth_phase{9.0 * 60.0};

  // Duration of apoptosis
  double duration_apoptosis{8.6 * 60.0};

  // Hypoxic threshold (numeric parameter \sigma_{H}, unit [1])
  double hypoxic_threshold{0.0538};

  // Probability parameters for the transition from Q to D modulated by the
  // nutrients N
  double threshold_Q_D_N{0.0538};
  double gamma_Q_D_N{0.0245 / 60.0};
  double alpha_Q_D_N{0.000408 / 60.0};
  double k_Q_D_N{50.0};

  // Probability parameters for the transition from Q to D modulated by the
  // drug DOX
  double zeta_Q_D_DOX{30};

  // Probability parameters for the transition from Q to D modulated by the
  // drug TRA
  double zeta_Q_D_TRA{30};

  // Probability parameters for the transition from G1 to SG2 modulated by the
  // nutrients N
  double threshold_Q_SG2_N{0.0538};
  double alpha_Q_SG2_N{0.0493 / 60.0};

  // Probability parameters for the transition from G1 to SG2 modulated by
  // the drug TRA
  // double threshold_Q_SG2_TRA{0.1};
  // double gamma_Q_SG2_TRA{1.0};
  double alpha_Q_SG2_TRA{5};
  // double k_Q_SG2_TRA{30.0};

  // Probability parameters for the transition from SG2 to SG2 modulated by the
  // drug DOX, e.g. probability to be trapped in SG2 due to drug DOX
  double threshold_SG2_SG2_DOX{0.1};
  double alpha_SG2_SG2_DOX{0.001};
  double k_SG2_SG2_DOX{30.0};

  // Probability parameters for the transition from SG2 to D modulated by the
  // drug DOX
  double threshold_SG2_D_DOX{0.1};
  double alpha_SG2_D_DOX{0.001};
  double k_SG2_D_DOX{30.0};

  // Probability parameters for the transition from H to D modulated by the
  // drugs DOX and TRA
  double base_rate_H_D{0.0001};
  double zeta_H_D_DOX{20};
  double zeta_H_D_TRA{20};

  // -----------------------------------------------------------------------
  // Agent-Continuum interaction parameters
  // -----------------------------------------------------------------------

  // Uptake rate of glucose by cells (unit [min^{-1}])
  double uptake_rate_glucose{0.0483 / 60.0};

  // Secretion rate tumor cells, i.e.how much VEGF is released by a tumor cell
  // per minute.
  double secretion_rate_vegf{0.03 / 60.0};

  // VEGF threshold for sprouting
  double vegf_threshold_sprouting{1e-3};

  // Nutrient supply by vessel (unit Nutrients / (Area * min)])
  double nutrient_supply_rate_vessel{0.0001};

  // VEGF consumption by vessel (unit Nutrients / (Area * min)])
  double vegf_consumption_rate_vessel{-0.0000};

  // Nutrient supply by vessel (unit Nutrients / (Area * min)])
  double dox_supply_rate_vessel{0.0};

  // Nutrient supply by vessel (unit Nutrients / (Area * min)])
  double tra_supply_rate_vessel{0.0};

  // Nutrient consumption by TumorCell (unit Nutrients / (min)])
  double nutrient_consumption_rate_tcell{-0.0001};

  // VEGF supply by TumorCell (unit Nutrients / (min)])
  double vegf_supply_rate_tcell{0.0001};

  // Nutrient supply by TumorCell (unit Nutrients / (min)])
  double dox_consumption_rate_tcell{-0.00000};

  // Nutrient supply by TumorCell (unit Nutrients / (min)])
  double tra_consumption_rate_tcell{-0.00000};

  // -----------------------------------------------------------------------
  // Forces
  // -----------------------------------------------------------------------

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

  // -----------------------------------------------------------------------
  // Continuum parameters
  // -----------------------------------------------------------------------

  // Resolution of the nutrients (glucose) diffusion grid. Gets forwarded to
  // DiffusionGrid constructor
  int diffusion_resolution_nutrients{50};

  // Initial value of the nutrients (glucose) concentration, uniform over grid
  // (unit [1])
  double initial_concentration_nutrients{0.5};

  // Diffusion coefficient for nutrients (glucose) (unit [\mu m / min])
  double diffusion_nutrients{50.0 / 60.0};

  // The decay constant of nutrients (glucose). It basically describes an
  // exponential decay for each point in the diffusion grid.
  double decay_rate_nutrients{0.00001};

  // Boundary condition for nutrients (glucose) diffusion grid. Gets forwarded
  // to DiffusionGrid constructor. This is the constant value that is set on the
  // boundary. Neumann and dirichlet boundary conditions are specified via the
  // experiment parameter.
  double boundary_condition_nutrients{0.0};

  // Resolution of the nutrients VEGF grid. Gets forwarded to
  // DiffusionGrid constructor
  int diffusion_resolution_vegf{50};

  // Initial value of the VEGF concentration, uniform over grid (unit [1])
  double initial_concentration_vegf{0.0};

  // Diffusion constant for VEGF
  double diffusion_vegf{40.0 / 60.0};

  // Decay constant for VEGF
  double decay_rate_vegf{0.0};

  // Boundary condition for VEGF diffusion grid. Gets forwarded to
  // DiffusionGrid constructor. This is the constant value that is set on the
  // boundary. Neumann and dirichlet boundary conditions are specified via the
  // experiment parameter.
  double boundary_condition_vegf{0.0};

  // Resolution of the nutrients TRA grid. Gets forwarded to
  // DiffusionGrid constructor
  int diffusion_resolution_tra{3};

  // Initial value of the TRA concentration, uniform over grid (unit [1])
  double initial_concentration_tra{0.0};

  // Diffusion constant for TRA
  double diffusion_tra{0.0};

  // Decay constant for TRA
  double decay_rate_tra{0.0};

  // Boundary condition for TRA diffusion grid. Gets forwarded to
  // DiffusionGrid constructor. This is the constant value that is set on the
  // boundary. Neumann and dirichlet boundary conditions are specified via the
  // experiment parameter.
  double boundary_condition_tra{0.0};

  // Resolution of the nutrients DOX grid. Gets forwarded to
  // DiffusionGrid constructor
  int diffusion_resolution_dox{3};

  // Initial value of the DOX concentration, uniform over grid (unit [1])
  double initial_concentration_dox{0.0};

  // Diffusion constant for DOX
  double diffusion_dox{0.0};

  // Decay constant for DOX
  double decay_rate_dox{0.0};

  // Boundary condition for DOX diffusion grid. Gets forwarded to
  // DiffusionGrid constructor. This is the constant value that is set on the
  // boundary. Neumann and dirichlet boundary conditions are specified via the
  // experiment parameter.
  double boundary_condition_dox{0.0};

  // -----------------------------------------------------------------------
  // Vessel parameters
  // -----------------------------------------------------------------------

  // Number of modes for random vessels
  int random_vessel_num_modes{10};

  // Exponent for random vessels
  double random_vessel_exponent{1.0};

  // Maximum deviation factor for random vessels (unit [1], 0 <= x <= 1).
  // Multiplies with vessel length to get the maximum deviation.
  double random_vessel_max_deviation{0.2};

  // Mean of the normal distribution for sinusoidal frequencies
  double random_vessel_frequency_mean{3.1415};

  // Standard deviation of the normal distribution for sinusoidal frequencies
  double random_vessel_frequency_std{4.0};

  // Length of vessels at initialization
  double default_vessel_length{5};

  // VEGF gradient threshold for apical growth
  double vegf_grad_threshold_apical_growth{1e-5};

  // Minimum distance to bifurcatoin or terminal end of vessel to allow
  // sprouting
  double min_dist_to_bifurcation{60.0};

  // Minimum distance to another tip cell to allow sprouting
  double min_dist_to_tip_cell{60.0};

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

  // Quotient threshold (stopping criterion) for apical growth
  double apical_growth_quotient_threshold{0.01};

  // -----------------------------------------------------------------------
  // Functions
  // -----------------------------------------------------------------------

  double GetSpheroidRadius() const {
    return std::pow(static_cast<double>(num_cells) * std::pow(cell_radius, 3) /
                        filled_volume,
                    1.0 / 3.0);
    ;
  }
};

}  // namespace bdm

#endif  // SIM_PARAM_H_
