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

#include "tumor_cell.h"
#include "core/environment/uniform_grid_environment.h"
#include "mechanical_forces.h"
#include "sim_param.h"
#include "transition_probabilities.h"

namespace bdm {

// Adds the first three components of y to x. The last component y[4] is
// ignored. Adds to x inplace.
inline void Add4dTo3dVector(Double3& x, const Double4& y) {
  x[0] += y[0];
  x[1] += y[1];
  x[2] += y[2];
}

void TumorCell::Initialize(const NewAgentEvent& event) {
  Base::Base::Initialize(event);  // Call the initialize of Agent not of Cell.

  if (event.GetUid() == CellDivisionEvent::kUid) {
    // The first few lines in the cell division are "shamelessly" copied /
    // adapted from biodynamo's cell.h.

    // First we get the pointers to the new and old TumorCell objects.
    const auto& cdevent = static_cast<const CellDivisionEvent&>(event);
    auto* mother_cell = bdm_static_cast<TumorCell*>(event.existing_agent);
    auto* daughter = this;

    // Inherit displacement scale factor
    daughter->displacement_scale_factor_ =
        mother_cell->displacement_scale_factor_;

    // This change is a littel problematic. Since we have to use diameter_ to
    // define the interaction range, we have to introduce additional variables
    // and change the very first line here. Hence, we cannot use
    // Cell::Initialize as a starting point but have to copy
    // double radius = mother_cell->GetDiameter() * 0.5;
    double radius = mother_cell->GetRadius();  // Change to cell.h

    // Define an axis for division (along which the cells will move). The angles
    // are coming from a random generator.
    double x_coord = std::cos(cdevent.theta) * std::sin(cdevent.phi);
    double y_coord = std::sin(cdevent.theta) * std::sin(cdevent.phi);
    double z_coord = std::cos(cdevent.phi);
    Double3 coords = {x_coord, y_coord, z_coord};
    double total_length_of_displacement = radius / displacement_scale_factor_;

    const auto& x_axis = mother_cell->kXAxis;
    const auto& y_axis = mother_cell->kYAxis;
    const auto& z_axis = mother_cell->kZAxis;

    Double3 axis_of_division =
        (coords.EntryWiseProduct(x_axis) + coords.EntryWiseProduct(y_axis) +
         coords.EntryWiseProduct(z_axis)) *
        total_length_of_displacement;

    // two equations for the center displacement :
    //  1) d2/d1= v2/v1 = volume_ratio (each sphere is shifted inver.
    //  proportionally to its volume)
    //  2) d1 + d2 = TOTAL_LENGTH_OF_DISPLACEMENT
    double d_2 = total_length_of_displacement / (cdevent.volume_ratio + 1);
    double d_1 = total_length_of_displacement - d_2;

    double mother_volume = mother_cell->GetVolume();
    double new_volume = mother_volume / (cdevent.volume_ratio + 1);
    daughter->SetVolume(mother_volume - new_volume);

    // position
    auto mother_pos = mother_cell->GetPosition();
    auto new_position = mother_pos + (axis_of_division * d_2);
    daughter->SetPosition(new_position);

    // E) This sphere becomes the 1st daughter
    // move these cells on opposite direction
    mother_pos -= axis_of_division * d_1;
    // update mother here and not in Update method to avoid recomputation
    mother_cell->SetPosition(mother_pos);
    mother_cell->SetVolume(new_volume);

    daughter->SetAdherence(mother_cell->GetAdherence());
    daughter->SetDensity(mother_cell->GetDensity());

    // SetVolume updates the diameter, we need to update the radius manually
    double new_radius_mother = std::cbrt(new_volume * 3 / (4 * Math::kPi));
    double new_radius_daughter =
        std::cbrt(daughter->GetVolume() * 3 / (4 * Math::kPi));
    mother_cell->SetRadius(new_radius_mother);
    daughter->SetRadius(new_radius_daughter);

    // Copy values of private member variables
    daughter->cell_state_ = mother_cell->cell_state_;
    daughter->t_last_state_transition_ = mother_cell->t_last_state_transition_;
    daughter->action_radius_factor_ = mother_cell->action_radius_factor_;
    daughter->growth_rate_ = mother_cell->growth_rate_;
    daughter->max_radius_ = mother_cell->max_radius_;
    daughter->nuclear_radius_ = mother_cell->nuclear_radius_;
    daughter->diffusion_substance_id_ = mother_cell->diffusion_substance_id_;

    // Update the action radi of the cells
    mother_cell->UpdateActionRadius();
    daughter->UpdateActionRadius();
  }
}

void TumorCell::SetRadii(double radius, double nuclear_radius,
                         double action_radius) {
  SetRadius(radius);
  SetNuclearRadius(nuclear_radius);
  SetActionRadius(action_radius);
}

inline void TumorCell::LimitDisplacementAtBoundary(Double3& displacement) {
  const auto* param = Simulation::GetActive()->GetParam();
  const auto* sparam = param->Get<SimParam>();
  const double min = sparam->lower_bound;
  const double max = sparam->upper_bound;
  const double r = this->radius_;
  const Double3 next_position = this->GetPosition() + displacement;
  for (int i = 0; i < 3; i++) {
    if (next_position[i] - r < min || next_position[i] + r > max) {
      displacement[i] = 0;
    }
  }
}

void TumorCell::ComputeApoptosisVolumeDecrease(double apoptosis_duration) {
  double lost_volume = 4.0 / 3.0 * Math::kPi *
                       (pow(nuclear_radius_, 3) - pow(radius_, 3)) /
                       apoptosis_duration;
  SetGrowthRate(lost_volume);
}

Double3 TumorCell::CalculateDisplacement(const InteractionForce* force,
                                         double squared_radius, double dt) {
  // 1. Get necessary objects for computation
  Double3 displacement{0, 0, 0};
  Double3 velocity{0, 0, 0};
  Double3 force_on_cell{0, 0, 0};

  auto* sim = Simulation::GetActive();
  auto* ctxt = sim->GetExecutionContext();
  const auto* param = sim->GetParam();
  const auto* sparam = param->Get<SimParam>();
  auto* env = dynamic_cast<UniformGridEnvironment*>(sim->GetEnvironment());

  // Set search radius manually to the one that was defined in hybrid_model.cc
  // This is necessary because we cannot use the one provided by the
  // MechanicalForcesOp, because it feeds pow(env->GetLargestAgent,2) to the
  // CalculateDisplacement Function.
  squared_radius = pow(static_cast<double>(env->GetBoxLength()), 2.0);

  // 2. Cast the force to our custom force
  const MechanicalInteractionForce* interaction_force =
      bdm_static_cast<const MechanicalInteractionForce*>(force);

  // 3. Iterate over all neighbours and compute forces onto this agent.
  auto calculate_neighbor_forces =
      L2F([&](Agent* neighbor, double squared_distance) {
        Double4 neighbor_force = interaction_force->Calculate(this, neighbor);
        Add4dTo3dVector(force_on_cell, neighbor_force);
      });
  ctxt->ForEachNeighbor(calculate_neighbor_forces, *this, squared_radius);

  // 4) The velocity of cell is the force divided by the viscosity (see Rocha
  // 2018 / Lima 2021). We limit the maximal speed with which a cell can move.
  velocity = force_on_cell / sparam->viscosity;
  double speed = velocity.Norm();
  if (speed > sparam->max_speed) {
    velocity *= sparam->max_speed / speed;
  }

  // 5) In Rocha 2018 / Lima 2021, they add additional forces from the domain
  // boundaries. Here, we simply cancel the movement in the direction of a
  // boundary if we are getting to close to the boundary.
  displacement = velocity * param->simulation_time_step;
  LimitDisplacementAtBoundary(displacement);
  return displacement;
}

void TumorCell::UpdateCellCycle() {
  auto* random = Simulation::GetActive()->GetRandom();
  const auto* param = Simulation::GetActive()->GetParam();
  const auto* sparam = param->Get<SimParam>();
  const double current_time =
      Simulation::GetActive()->GetScheduler()->GetSimulatedTime();
  const double time_in_state = current_time - t_last_state_transition_;
  const double duplication_time =
      sparam->duration_cell_cycle - sparam->duration_growth_phase;
  // Quiescent states stochastically transition into the states proliferative
  // or dead with certain probabilities.
  if (cell_state_ == CellState::kQuiescent) {
    double sigma{0.0};
    auto* diffusion_grid =
        Simulation::GetActive()->GetResourceManager()->GetDiffusionGrid(
            diffusion_substance_id_);
    sigma = diffusion_grid->GetConcentration(GetPosition());
    double probability_death =
        ComputeProbabilityDeath(sigma, param->simulation_time_step, sparam);
    double probability_prolif = ComputeProbabilityProliferative(
        sigma, param->simulation_time_step, sparam);
    double decision_variable = random->Uniform();
    if (decision_variable < probability_prolif) {
      SetCellState(CellState::kProliferativeSG2);
      t_last_state_transition_ = current_time;
    } else if (decision_variable < probability_prolif + probability_death) {
      SetCellState(CellState::kDead);
      t_last_state_transition_ = current_time;
      // decrease of cell volume per step due to apoptosis
      ComputeApoptosisVolumeDecrease(sparam->duration_apoptosis);
    } else {
      ;
    }
    // Cells in the state kProliferativeSG2 don't do anything until the have
    // finished duplicating their DNA (etc.) after \tau_p-\tau_G1. Once that
    // has happened, the cell divides and enters into the G1 growth phase.
  } else if (cell_state_ == CellState::kProliferativeSG2 &&
             time_in_state > duplication_time) {
    t_last_state_transition_ = current_time;
    SetCellState(CellState::kProliferativeG1);
    Divide();
    // Cells in the cell states kProliferativeG1 increase their volume linearly
    // and after a certain time \tau_G1 they stop and enter the quiescent state.
  } else if (cell_state_ == CellState::kProliferativeG1) {
    // ToDo increase volume linearly;
    ChangeVolume(growth_rate_);
    if (time_in_state > sparam->duration_growth_phase) {
      SetCellState(CellState::kQuiescent);
      t_last_state_transition_ = current_time;
    }
  } else if (cell_state_ == CellState::kDead) {
    // decrease volume with previously computed negative growth rate
    // (ComputeApoptosisVolumeDecrease)
    if (!sparam->keep_dead_cells) {
      ChangeVolume(growth_rate_);
    }
  } else {  // In any other case do nothing.
    ;       // do nothing
  }
}

void TumorCell::ChangeVolume(double speed) {
  // scaling for integration step
  auto* param = Simulation::GetActive()->GetParam();
  double delta = speed * param->simulation_time_step;
  double volume = GetVolume();
  volume += delta;
  // Compute the new redius of the cell from volume
  double radius = std::cbrt(volume * 3 / (4 * Math::kPi));
  if (radius > radius_) {
    Base::SetPropagateStaticness();  // copied from cell implementation
  }
  if (volume < 5.2359877E-7) {
    Log::Fatal("TumorCell::ChangeVolume", "Cell volume is getting too small.");
  }
  // Set the new radius and update Diameter+ActionRadius.
  if (radius < max_radius_) {
    SetRadius(radius);
    UpdateActionRadius();
    // Set the volume to the new value. ToDo / WARNING: this function is not
    // Available in bdm master, it's a dirty work around. Find better way to do
    // that.
    SetVolume(volume);
  } else {
    Log::Warning("TumorCell::ChangeVolume",
                 "Cell has reached maximal possible size");
  }
}

void ProgressInCellCycle::Run(Agent* agent) {
  auto* tumor_cell = bdm_static_cast<TumorCell*>(agent);
  if (tumor_cell != nullptr) {
    tumor_cell->UpdateCellCycle();
  }
}

void Death::Run(Agent* agent) {
  auto* tumor_cell = bdm_static_cast<TumorCell*>(agent);
  if (tumor_cell->GetCellState() == CellState::kDead &&
      tumor_cell->GetRadius() < tumor_cell->GetNuclearRadius()) {
    tumor_cell->RemoveFromSimulation();
  }
}

}  // namespace bdm
