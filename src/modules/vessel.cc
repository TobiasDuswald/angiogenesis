// //
// -----------------------------------------------------------------------------
// //
// // Copyright (C) 2022 CERN, TUM, and UT Austin. All Rights Reserved.
// //
// // Licensed under the Apache License, Version 2.0 (the "License");
// // you may not use this file except in compliance with the License.
// //
// // See the LICENSE file distributed with this work for details.
// // See the NOTICE file distributed with this work for additional information
// // regarding copyright ownership.
// //
// //
// -----------------------------------------------------------------------------

#include "vessel.h"
#include "modules/tumor_cell.h"

namespace bdm {

void SproutingAngiogenesis::Initialize(const NewAgentEvent& event) {
  Base::Initialize(event);
  can_branch_ = false;
}

void SproutingAngiogenesis::Run(Agent* agent) {
  auto* sim = Simulation::GetActive();
  auto* random = sim->GetRandom();
  auto* rm = sim->GetResourceManager();
  // const auto* const sparam = sim->GetParam()->Get<SimParam>();

  // First time the behaviour is executed we get and remember a pointer to the
  // relevant diffusion grid
  if (!init_) {
    dg_guide_ = rm->GetDiffusionGrid(Substances::kVEGF);
    init_ = true;
  }

  // Downcast agent to Vessel
  auto* dendrite = dynamic_cast<Vessel*>(agent);

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
  double min_distance_to_bifurcation = 60.0;
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

    // Compute sprouting direction on cone around gradient
    double phi = random->Uniform(2 * Math::kPi);
    double theta = random->Uniform(0.25, 0.80);
    Double3 sprouting_direction = VectorOnConeAroundAxis(gradient, phi, theta);

    // Branch vessel
    // auto* dendrite2 =
    dendrite->Branch(dendrite->GetDiameter(), sprouting_direction,
                     dendrite->GetActualLength() / 2);
    // auto* casted_vessel = dynamic_cast<Vessel*>(dendrite2);
    // casted_vessel->AddBehavior(new NutrientSupply("Nutrients", 1.0));
  }
}

void ApicalGrowth::Initialize(const NewAgentEvent& event) {
  Base::Initialize(event);
  can_branch_ = false;
}

void ApicalGrowth::Run(Agent* agent) {
  auto* sim = Simulation::GetActive();
  auto* random = sim->GetRandom();
  auto* rm = sim->GetResourceManager();
  // const auto* const sparam = sim->GetParam()->Get<SimParam>();

  if (!init_) {
    dg_guide_ = rm->GetDiffusionGrid(Substances::kVEGF);
    init_ = true;
  }

  // Cast Agent to dendrite
  auto* dendrite = dynamic_cast<Vessel*>(agent);

  /// 1. Check if element can grow
  if (!dendrite || !dendrite->IsTerminal() || !dendrite->CanGrow()) {
    return;
  }

  /// 2. Get gradient and check if magnitude of gradient is above a
  ///    threshold (0 to begin with).
  double threshold_gradient{0.0001};
  Double3 gradient;
  dg_guide_->GetGradient(dendrite->GetPosition(), &gradient, false);
  if (gradient.Norm() < threshold_gradient) {
    return;
  }

  /// 3. If vessel is close to an Tumor cell, we interrupt the growth as well
  double distance_for_growth_stop = 40;
  CountNeighborsFunctor<TumorCell*> cnf(distance_for_growth_stop);
  auto* ctxt = Simulation::GetActive()->GetExecutionContext();
  ctxt->ForEachNeighbor(cnf, *dendrite, std::pow(distance_for_growth_stop, 2));
  if (cnf.GetNumNeighbors() != 0) {
    dendrite->ProhibitGrowth();
    return;
  }

  /// 4. Extend into gradient direction with some random disturbance and
  ///    memory.
  double weight_random{0.2};
  double weight_old{0.4};
  double weight_gradient{0.4};
  double growth_speed{1.0};
  auto random_direction = random->template UniformArray<3>(-1, 1);
  auto old_direction = dendrite->GetSpringAxis();

  Double3 new_direction = old_direction * weight_old +
                          random_direction * weight_random +
                          gradient.GetNormalizedArray() * weight_gradient;

  dendrite->ElongateTerminalEnd(growth_speed, new_direction);
}

void NutrientSupply::Initialize(const NewAgentEvent& event) {
  Base::Initialize(event);
  auto* other = dynamic_cast<NutrientSupply*>(event.existing_behavior);
  substance_ = other->substance_;
  dgrid_ = other->dgrid_;
  quantity_ = other->quantity_;
}

void NutrientSupply::Run(Agent* agent) {
  auto* vessel = dynamic_cast<Vessel*>(agent);
  // If the behaviour is assigned to a vessel and it is not part of the initial
  // vascular network, we do supply the nutrients.
  if (vessel && vessel->CanGrow()) {
    // This is very simple way to supply nutrients. In the current setup,
    // bifurcations add more nutrients than other elements. This should probably
    // be changed. Also note that the behavior makes sense if the vessel is
    // in the same scale as the diffusion grid, or smaller.

    // Weights for end, middle and start of the vessel, respectively.
    std::array<double, 3> weights = {0.25, 0.5, 0.25};

    // Collect positions (max 3) in vector. Necessary because theoretically, a
    // vessel can have a soma as mother compartment.
    std::vector<Double3> positions;

    // Collect positions of the vessel.
    const auto& end = vessel->GetMassLocation();
    const auto& middle = vessel->GetPosition();
    positions.push_back(end);
    positions.push_back(middle);
    auto* mother_ptr = dynamic_cast<Vessel*>(vessel->GetMother().Get());
    if (mother_ptr != nullptr) {
      const auto& start = mother_ptr->GetMassLocation();
      positions.push_back(start);
    }
    for (size_t i = 0; i < positions.size(); i++) {
      // Change concentration according to the weight.
      dgrid_->ChangeConcentrationBy(positions[i], quantity_ * weights[i]);
    }
  }
}
}  // namespace bdm
