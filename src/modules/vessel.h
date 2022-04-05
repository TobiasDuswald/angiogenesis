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

namespace bdm {

class Vessel : public NeuriteElement {
  BDM_AGENT_HEADER(Vessel, NeuriteElement, 1);

 public:
  Vessel() : can_grow_(true){};

  void AllowGrowth() { can_grow_ = true; }
  void ProhibitGrowth() { can_grow_ = false; }
  bool CanGrow() { return can_grow_; }

  void RunDiscretization() override {
    if (!CanGrow()) {
      // For vessel agents that are part of the initial vasculature, we do not
      // execute the Discretization() function. The discretization generates new
      // vessels, which are allowed to grow and therefore also secrete
      // nutrients. We do not want our initial vasculature to supply nutrients.
      return;
    }
    Base::RunDiscretization();
  }

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
  /// Boolean to execute certain calls only upon initialization
  bool init_ = false;
  bool can_branch_ = true;
  /// DiffusionGrid for guiding substance
  DiffusionGrid* dg_guide_ = nullptr;
};

/// Behaviour to grow vessels towards higher VEGF concentrations.
class ApicalGrowth : public Behavior {
  BDM_BEHAVIOR_HEADER(ApicalGrowth, Behavior, 1);
  ApicalGrowth() { AlwaysCopyToNew(); }
  virtual ~ApicalGrowth() {}

  void Initialize(const NewAgentEvent& event) override;

  void Run(Agent* agent) override;

 private:
  bool init_ = false;
  bool can_branch_ = true;
  DiffusionGrid* dg_guide_ = nullptr;
};

/// Supply nutrients to surrounding tissues.
class NutrientSupply : public Behavior {
  BDM_BEHAVIOR_HEADER(NutrientSupply, Behavior, 1);

 public:
  NutrientSupply() { AlwaysCopyToNew(); };
  explicit NutrientSupply(const std::string& substance, double quantity = 1)
      : substance_(substance), quantity_(quantity) {
    dgrid_ = Simulation::GetActive()->GetResourceManager()->GetDiffusionGrid(
        substance);
    AlwaysCopyToNew();
  }

  explicit NutrientSupply(DiffusionGrid* dgrid, double quantity = 1)
      : dgrid_(dgrid), quantity_(quantity) {
    substance_ = dgrid->GetSubstanceName();
  }

  virtual ~NutrientSupply() = default;

  void Initialize(const NewAgentEvent& event) override;

  void Run(Agent* agent) override;

 protected:
  std::string substance_;
  DiffusionGrid* dgrid_ = nullptr;
  double quantity_ = 1;
};

}  // namespace bdm

#endif  // VESSEL_H_
