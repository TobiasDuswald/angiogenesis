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

  void SetRadii(double radius, double nuclear_radius, double action_radius);
};

// Behaviour that allows cells to transition between states and progress in the
// cell cycle. It also contains growth, shrinking, and the cell division.
struct ProgressInCellCycle : public Behavior {
  BDM_BEHAVIOR_HEADER(ProgressInCellCycle, Behavior, 1);

  ProgressInCellCycle() { AlwaysCopyToNew(); }

  void Run(Agent* agent) override;
};

// A simple behaviour that turns a TumorCell hypoxic or quiescent depending on
// the concentration of substance_id (threshold behavior).
struct UpdateHypoxic : public Behavior {
  BDM_BEHAVIOR_HEADER(UpdateHypoxic, Behavior, 1);

 public:
  UpdateHypoxic() = default;
  explicit UpdateHypoxic(int substance_id) : substance_id_(substance_id) {
    AlwaysCopyToNew();
  }

  void Run(Agent* agent) override;

 protected:
  int substance_id_ = 0;
};

/// Secrete a substance only if it's cell state is hypoxic.
// ToDo: Implement more consciese by inheriting from Secretion.
class HypoxicSecretion : public Behavior {
  BDM_BEHAVIOR_HEADER(HypoxicSecretion, Behavior, 1);

 public:
  HypoxicSecretion() = default;
  explicit HypoxicSecretion(const std::string& substance, double quantity = 1)
      : substance_(substance), quantity_(quantity) {
    dgrid_ = Simulation::GetActive()->GetResourceManager()->GetDiffusionGrid(
        substance);
  }

  explicit HypoxicSecretion(DiffusionGrid* dgrid, double quantity = 1)
      : dgrid_(dgrid), quantity_(quantity) {
    substance_ = dgrid->GetContinuumName();
  }

  virtual ~HypoxicSecretion() = default;

  void Initialize(const NewAgentEvent& event) override;

  void Run(Agent* agent) override;

 protected:
  std::string substance_;
  DiffusionGrid* dgrid_ = nullptr;
  double quantity_ = 1;
};

/// Supply nutrients to surrounding tissues. The vessel is discretized along
/// it's center axis into N points. The number N is computed automatically such
/// that we have roughly 2 points per voxel. When using the behavior, we specify
/// a quantity - note that this quantity is weighted with the vessel-agent's
/// surface and corrected by a term that avoids overshooting the maximum
/// concentration (logistic growth).
class PointContinuumInteraction : public Behavior {
  BDM_BEHAVIOR_HEADER(PointContinuumInteraction, Behavior, 1);

 public:
  PointContinuumInteraction() { AlwaysCopyToNew(); };
  explicit PointContinuumInteraction(double rate_nutrients, double rate_vegf,
                                     double rate_dox, double rate_tra)
      : interaction_rate_({rate_nutrients, rate_vegf, rate_dox, rate_tra}) {
    // Always copy the behavior to a new agent.
    AlwaysCopyToNew();
  }

  virtual ~PointContinuumInteraction() = default;

  void Initialize(const NewAgentEvent& event) override;

  void Run(Agent* agent) override;

 private:
  std::array<double, 4> interaction_rate_ = {0, 0, 0, 0};
  bool init_ = false;
};

}  // namespace bdm

#endif  // TUMOR_CELL_H_
