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
  Vessel() : can_grow_(true){};

  void AllowGrowth() { can_grow_ = true; }
  void ProhibitGrowth() { can_grow_ = false; }
  bool CanGrow() { return can_grow_; }

  /// Method called by default discretization operation
  void RunDiscretization() override;

  /// Returns the surface area of the cylinder
  double GetSurfaceArea() const;

 protected:
  /// Parameter to decide if a vessel compartment can grow towards a higher
  /// VEGF concentration (used to fix initial vessel configuration)
  bool can_grow_;
};

/// Behaviour to create a new bifurcation if external VEGF concentarion is
/// surpassing a threshold.
class SproutingAngiogenesis : public Behavior {
  BDM_BEHAVIOR_HEADER(SproutingAngiogenesis, Behavior, 1);
  SproutingAngiogenesis() { AlwaysCopyToNew(); }
  virtual ~SproutingAngiogenesis() {}

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
  virtual ~ApicalGrowth() {}

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
class NutrientSupply : public Behavior {
  BDM_BEHAVIOR_HEADER(NutrientSupply, Behavior, 1);

 public:
  NutrientSupply() { AlwaysCopyToNew(); };
  explicit NutrientSupply(const std::string& substance, double quantity = 1,
                          bool use_for_consumption = false)
      : quantity_(quantity), use_for_consumption_(use_for_consumption) {
    // Register Pointer to substance
    dgrid_ = Simulation::GetActive()->GetResourceManager()->GetDiffusionGrid(
        substance);
    // Always copy the behavior to a new agent.
    AlwaysCopyToNew();
  }

  virtual ~NutrientSupply() = default;

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

  DiffusionGrid* dgrid_ = nullptr;
  size_t n_sample_points_ = 3;
  double quantity_ = 1;

  // Track if we have already figured out the
  bool init_ = false;
  // Do we recycle this behavior to consume nutrients or VEGF?
  bool use_for_consumption_ = false;
};

}  // namespace bdm

#endif  // VESSEL_H_
