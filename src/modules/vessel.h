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

namespace bdm {

class Vessel : public NeuriteElement {
  BDM_AGENT_HEADER(Vessel, NeuriteElement, 1);

 public:
  Vessel() : can_grow_(true){};

  void AllowGrowth() { can_grow_ = true; }
  void ProhibitGrowth() { can_grow_ = false; }
  bool CanGrow() { return can_grow_; }

 protected:
  /// Parameter to decide if a vessel compartment can grow towards a higher VEGF
  /// concentration (used to fix initial vessel configuration)
  bool can_grow_;
};

/// Behaviour to create a new bifurcation if external VEGF concentarion is
/// surpassing a threshold.
class SproutingAngiogenesis : public Behavior {
  BDM_BEHAVIOR_HEADER(SproutingAngiogenesis, Behavior, 1);
  SproutingAngiogenesis() { AlwaysCopyToNew(); }
  virtual ~SproutingAngiogenesis() {}

  void Initialize(const NewAgentEvent& event) override {
    Base::Initialize(event);
    can_branch_ = false;
  }

  void Run(Agent* agent) override {
    auto* sim = Simulation::GetActive();
    auto* random = sim->GetRandom();
    auto* rm = sim->GetResourceManager();
    const auto* const sparam = sim->GetParam()->Get<SimParam>();

    // First time the behaviour is executed we get and remember a pointer to the
    // relevant diffusion grid
    if (!init_) {
      dg_guide_ = rm->GetDiffusionGrid(Substances::kVEGF);
      init_ = true;
    }

    // Downcast agent to NeuriteElement
    auto* dendrite = dynamic_cast<NeuriteElement*>(agent);

    /// Check if vessel can branch
    if (dendrite->GetDaughterLeft() == nullptr ||
        dendrite->GetDaughterRight() != nullptr) {
      return;
    }

    /// 1. Check if the concentration of VEGF is above a certain threshold
    double my_vegf_threshold = 1e-3;
    double my_vegf_concentration =
        dg_guide_->GetConcentration(dendrite->GetPosition());
    if (my_vegf_concentration < my_vegf_threshold) {
      return;
    }

    /// 2. Check if the is eligible for branching, e.g. minium distance to other
    ///    bifurcation points is given.
    double min_distance_to_bifurcation = 1.0;
    double distance{0.0};
    // 2.1 Walk down the tree ("daughter" direction)
    AgentPointer<NeuriteElement> daughter = dendrite->GetDaughterLeft();
    while (true) {
      distance += daughter->GetActualLength();
      if (distance >= min_distance_to_bifurcation) {
        // If we keep the minumum distance to the next bifurcation, we are happy
        // end exit the loop.
        break;
      } else if (daughter->GetDaughterLeft() && !daughter->GetDaughterRight()) {
        // If we can follow the vessel further and we do not reach a
        // bifurcation, we take the next vessel element.
        daughter = daughter->GetDaughterLeft();
      } else {
        // In case we hit the end of a vessel or a bifurcation without having
        // achieved a minimum distance, we end the behaviour here.
        return;
      }
    }
    // 2.2 Walk up the tree ("mother" direction)
    distance = 0.0;
    AgentPointer<neuroscience::NeuronOrNeurite> mother = dendrite->GetMother();
    while (true) {
      auto* mother_ptr = mother.Get();
      if (mother_ptr == nullptr) {
        Log::Fatal("SproutingAngiogenesis::Run", "mother is null.");
      } else if (mother->IsNeuriteElement()) {
        auto* mother_neurite_ptr = static_cast<NeuriteElement*>(mother_ptr);
        if (mother_neurite_ptr->GetDaughterRight()) {
          // If we reach a bifurcation / branching point before reaching
          // min_distance_to_bifurcation we stop.
          return;
        }
        distance += mother_neurite_ptr->GetActualLength();
        mother = mother_neurite_ptr->GetMother();
      } else {
        // Mother is soma, i.e. end of vessel. Then we have not reached the
        // necessary distance and exit
        return;
      }
      if (distance >= min_distance_to_bifurcation) {
        // If we keep the minumum distance to the next bifurcation, we are happy
        // end exit the loop.
        break;
      }
    }

    /// 3. If both criteria are satisfied (i.e. we reach this point), create a
    ///    sprout with a certain probability growing towards the gradient of
    ///    VEGF (minimum angle).
    double sprouting_probability{0.001};
    if (random->Uniform() < sprouting_probability) {
      // Get gradient
      Double3 gradient;
      dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);

      // Compute sprouting direction
      Double3 sprouting_direction = gradient;

      // Branch vessel
      auto* dendrite_2 =
          dendrite->Branch(dendrite->GetDiameter(), sprouting_direction,
                           dendrite->GetActualLength() / 2);
    }
  }

 private:
  /// Boolean to execute certain calls only upon initialization
  bool init_ = false;
  bool can_branch_ = true;
  /// DiffusionGrid for guiding substance
  DiffusionGrid* dg_guide_ = nullptr;
};

}  // namespace bdm

#endif  // VESSEL_H_
