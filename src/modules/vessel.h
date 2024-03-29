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

#ifndef VESSEL_H_
#define VESSEL_H_

#include "biodynamo.h"
#include "neuroscience/neuroscience.h"
#include "sim_param.h"
#include "util/neighbor_counter.h"
#include "util/vector_operations.h"

#include <algorithm>

namespace bdm {

class Vessel : public NeuriteElement {
  BDM_AGENT_HEADER(Vessel, NeuriteElement, 1);

 public:
  Vessel() = default;

  void AllowGrowth() { can_grow_ = true; }
  void ProhibitGrowth() { can_grow_ = false; }
  bool CanGrow() const { return can_grow_; }

  /// Method called by default discretization operation
  void RunDiscretization() override;

  /// Deactivate displacement calculation. Our force does not consider forces
  /// on the vessel itself, but the displacement calculation triggers a
  /// ForEachNeighbor call, which is not necessary. This effectively saves us
  /// one ForEachNeighbor call per vessel agent.
  Real3 CalculateDisplacement(const InteractionForce* force,
                              real_t squared_radius, real_t dt) override {
    return Real3({0, 0, 0});
  }

  /// Returns if the vessel is a tip cell
  bool IsTipCell() const;

  /// Returns if the vessel is a stalk cell
  bool IsStalkCell() const;

  /// Returns the surface area of the cylinder
  double GetSurfaceArea() const;

 private:
  /// Parameter to decide if a vessel compartment can grow towards a higher
  /// VEGF concentration (used to fix initial vessel configuration)
  bool can_grow_{true};

  // Split the vessel into two parts at the given position
  Vessel* SplitVessel(real_t distal_portion);
};

/// Behaviour to create a new bifurcation if external VEGF concentarion is
/// surpassing a threshold.
class SproutingAngiogenesis : public Behavior {
  BDM_BEHAVIOR_HEADER(SproutingAngiogenesis, Behavior, 1);
  SproutingAngiogenesis() { AlwaysCopyToNew(); }
  ~SproutingAngiogenesis() override = default;

  void Initialize(const NewAgentEvent& event) override;

  void Run(Agent* agent) override;

 private:
  /// DiffusionGrid for guiding substance
  DiffusionGrid* dg_guide_ = nullptr;
  /// Boolean to execute certain calls only upon initialization
  bool init_ = false;
  bool can_branch_ = true;
};

/// Behaviour to grow vessels towards higher VEGF concentrations.
class ApicalGrowth : public Behavior {
  BDM_BEHAVIOR_HEADER(ApicalGrowth, Behavior, 1);
  ApicalGrowth() { AlwaysCopyToNew(); }
  ~ApicalGrowth() override = default;

  void Initialize(const NewAgentEvent& event) override;

  void Run(Agent* agent) override;

 private:
  DiffusionGrid* dg_guide_ = nullptr;
  bool init_ = false;
  bool can_branch_ = true;
};

/// Supply nutrients to surrounding tissues. The vessel is discretized along
/// it's center axis into N points. The number N is computed automatically such
/// that we have roughly 2 points per voxel. When using the behavior, we specify
/// a quantity - note that this quantity is weighted with the vessel-agent's
/// surface and corrected by a term that avoids overshooting the maximum
/// concentration (logistic growth).
class LineContinuumInteraction : public Behavior {
  BDM_BEHAVIOR_HEADER(LineContinuumInteraction, Behavior, 1);

 public:
  LineContinuumInteraction() { AlwaysCopyToNew(); };
  explicit LineContinuumInteraction(double rate_nutrients, double rate_vegf,
                                    double rate_dox, double rate_tra)
      : interaction_rate_({rate_nutrients, rate_vegf, rate_dox, rate_tra}) {
    // Always copy the behavior to a new agent.
    AlwaysCopyToNew();
  }

  ~LineContinuumInteraction() override = default;

  void Initialize(const NewAgentEvent& event) override;

  void Run(Agent* agent) override;

 private:
  /// Compute weights for the sampling points. Weights are designed to add up
  /// to 1.0 and that the boundary points are weighted half as much as the
  /// interior points.
  void ComputeWeights();

  /// Compute the sample points for the source term along the center line.
  void ComputeSamplePoints(const Double3& start, const Double3& end);

  // Note, it's not very elegant to have a vector of pointers to doubles here.
  // However, it is the easiest way to implement the weights. The storage is
  // repetitive.
  std::vector<double> sample_weights_;
  std::vector<Double3> sample_points_;

  std::array<double, 4> interaction_rate_ = {0, 0, 0, 0};
  size_t n_sample_points_ = 0;
  double smallest_voxel_size_ = 10;

  // Track if we have already figured out the
  bool init_ = false;
};

}  // namespace bdm

#endif  // VESSEL_H_
