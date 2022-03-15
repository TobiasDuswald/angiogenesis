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

#ifndef TUMOR_CELL_H_
#define TUMOR_CELL_H_

#include "biodynamo.h"

namespace bdm {

// To label the states of a cell, we use the enum CellState.
enum CellState {
  kQuiescent,
  kProliferativeSG2,
  kProliferativeG1,
  kHypoxic,
  kDead
};

// Our class extends the Cell object
class TumorCell : public Cell {
  BDM_AGENT_HEADER(TumorCell, Cell, 1);

 private:
  // Cell state representing, for instance, quiescent, proliferative, or dead.
  int cell_state_;
  // The cell needs to know to which diffusion grid it is coupled.
  int diffusion_substance_id_;
  // Time of the last state transition.
  double t_last_state_transition_;
  // Radius of the cell.
  double radius_;
  // Radius of the nucleus.
  double nuclear_radius_;
  // Action radius of the cell.
  double action_radius_;
  // Action radius factor: action_radius_ = action_radius_factor_ * radius_
  double action_radius_factor_;
  // Growth rate of this specific cell.
  double growth_rate_;
  // Max radius that can be achieved for this cell.
  double max_radius_;
  // Factor to modify the displacement.
  double displacement_scale_factor_;

 public:
  TumorCell() {}
  explicit TumorCell(const Double3& position, int cell_state)
      : Base(position), cell_state_{cell_state} {
    action_radius_factor_ = 1.214;
    growth_rate_ = 0;
    max_radius_ = 20.0;
    displacement_scale_factor_ = 4.0;
    t_last_state_transition_ = 0.0;
    diffusion_substance_id_ = 0;
  }
  virtual ~TumorCell() {}

  // If TumorCell divides, the daughter has to initialize its attributes. This
  // member is called in that case.
  void Initialize(const NewAgentEvent& event) override;

  // Calculate the displacement of the TumorCell.
  Double3 CalculateDisplacement(const InteractionForce* force,
                                double squared_radius, double dt) override;

  // This function is called by the function CalculateDisplacement and prevents
  // cells leaving the defined simulation boundaries
  void LimitDisplacementAtBoundary(Double3& displacement);

  // This member function implements the stochastic and deterministic cell state
  // transitions as outlined in Lima 2021.
  void UpdateCellCycle();

  // This function is defined in class Cell but we overwrite it here to obtain
  // the correct behavior and increase all radii correctly.
  void ChangeVolume(double speed);

  // When entering the Apoptosis phase, cells start decreasing the volume. This
  // function changes the variable growth_rate_ to a negative value such that
  // the cell shrinks to the size of the nucleus in the defined time
  // apoptosis_duration.
  void ComputeApoptosisVolumeDecrease(double apoptosis_duration);

  //////////////////////////////////////////////////////////////////////////////
  // Getter and Setter functions
  //////////////////////////////////////////////////////////////////////////////
  void SetCellState(int cell_state) { cell_state_ = cell_state; }
  int GetCellState() { return cell_state_; };

  void SetRadius(double radius) {
    radius_ = radius;
    SetVolume(4.0 / 3.0 * Math::kPi * pow(radius, 3));
  }
  double GetRadius() const { return radius_; }

  void SetNuclearRadius(double nuclear_radius) {
    nuclear_radius_ = nuclear_radius;
  }
  double GetNuclearRadius() const { return nuclear_radius_; }

  void SetActionRadius(double action_radius) { action_radius_ = action_radius; }
  double GetActionRadius() const { return action_radius_; }

  void UpdateActionRadius() {
    SetActionRadius(action_radius_factor_ * radius_);
  }

  void SetActionRadiusFactor(double action_radius_factor) {
    action_radius_factor_ = action_radius_factor;
  }
  double GetActionRadiusFactor() { return action_radius_factor_; }

  void SetGrowthRate(double growth_rate) { growth_rate_ = growth_rate; }
  double GetGrowthRate() { return growth_rate_; }

  void SetMaxRadius(double max_radius) { max_radius_ = max_radius; }
  double GetMaxRadius() { return max_radius_; }

  void SetDisplacementScaleFactor(double displacement_scale_factor) {
    displacement_scale_factor_ = displacement_scale_factor;
  }
  double GetDisplacementScaleFactor() { return displacement_scale_factor_; }

  void SetSubstanceID(int id) { diffusion_substance_id_ = id; }

  void SetRadii(double radius, double nuclear_radius, double action_radius);
};

// Behaviour that allows cells to transition between states and progress in the
// cell cycle. It also contains growth, shrinking, and the cell division.
struct ProgressInCellCycle : public Behavior {
  BDM_BEHAVIOR_HEADER(ProgressInCellCycle, Behavior, 1);

  ProgressInCellCycle() { AlwaysCopyToNew(); }

  void Run(Agent* agent) override;
};

// Behaviour that allows cells to die, i.e. be removed from the simulation if
// the cell's volume has decreased sufficiently and it's in the dead state.
struct Death : public Behavior {
  BDM_BEHAVIOR_HEADER(Death, Behavior, 1);

  Death() { AlwaysCopyToNew(); }

  void Run(Agent* agent) override;
};

}  // namespace bdm

#endif  // TUMOR_CELL_H_
