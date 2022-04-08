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
  const auto* const sparam = sim->GetParam()->Get<SimParam>();

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
  double my_vegf_concentration =
      dg_guide_->GetConcentration(dendrite->GetPosition());
  if (my_vegf_concentration < sparam->vegf_threshold_sprouting) {
    return;
  }

  /// 2. Check if the is eligible for branching, e.g. minium distance to other
  ///    bifurcation points is given.
  double min_distance_to_bifurcation = sparam->min_dist_to_bifurcation;
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
  if (random->Uniform() < sparam->sprouting_probability) {
    // Get gradient
    Double3 gradient;
    dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);

    // Compute sprouting direction on cone around gradient
    double phi = random->Uniform(2 * Math::kPi);
    double theta = random->Uniform(0.25, 0.80);
    Double3 sprouting_direction = VectorOnConeAroundAxis(gradient, phi, theta);

    // Branch vessel
    dendrite->Branch(dendrite->GetDiameter(), sprouting_direction,
                     dendrite->GetActualLength() / 2);
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
  const auto* const sparam = sim->GetParam()->Get<SimParam>();

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
  Double3 gradient;
  dg_guide_->GetGradient(dendrite->GetPosition(), &gradient, false);
  if (gradient.Norm() < sparam->vegf_grad_threshold_apical_growth) {
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
  double weight_random = sparam->apical_growth_random_weight;
  double weight_old = sparam->apical_growth_old_weight;
  double weight_gradient = sparam->apical_growth_gradient_weight;
  double growth_speed = sparam->apical_growth_speed;
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

    bool skip_first_weight = false;
    auto* mother_ptr = dynamic_cast<Vessel*>(vessel->GetMother().Get());
    if (mother_ptr &&
        (vessel->GetAgentPtr<Vessel>() == mother_ptr->GetDaughterRight())) {
      skip_first_weight = true;
    }

    // Get start and and point
    const auto end = vessel->GetMassLocation();
    Double3 start;
    if (mother_ptr) {
      start = mother_ptr->GetMassLocation();
    } else {
      // If mother is soma, start with middle of vessel
      start = vessel->GetPosition();
    }

    if (!init_) {
      // In the first step, we define how granular the supply is. We do this by
      // demanding that the distance between two sampling points is roughly half
      // the box length of the discretization. Note that we compute the
      // the number of sample points here and not in the constructor because
      // the BoxLength is not initialized in the DiffusionGrid constructor.
      init_ = true;
      double delta_x = dgrid_->GetBoxLength();
      double distance = (end - start).Norm();
      n_sample_points_ =
          std::max(3, static_cast<int>(std::ceil(2 * distance / delta_x + 1)));
      ComputeWeights();
    }

    // Compute the sample points
    ComputeSamplePoints(start, end);

    // Add nutrients to the continuum.
    for (size_t i = 0; i < n_sample_points_; i++) {
      // We skip the first weight if we are in the right daughter vessel.
      if (skip_first_weight && i == 0) {
        continue;
      }
      // Change concentration according to the weight.
      double delta_concentration =
          quantity_ * sample_weights_[i] * vessel->GetSurfaceArea();
      dgrid_->ChangeConcentrationBy(sample_points_[i], delta_concentration);
    }
  }
}

void NutrientSupply::ComputeWeights() {
  sample_weights_.resize(n_sample_points_);
  for (size_t i = 0; i < n_sample_points_; i++) {
    sample_weights_[i] = 1 / static_cast<double>(n_sample_points_ - 1);
  }
  sample_weights_[0] *= 0.5;
  sample_weights_[n_sample_points_ - 1] *= 0.5;
}

void NutrientSupply::ComputeSamplePoints(const Double3& start,
                                         const Double3& end) {
  // Computes the sample points distributed along the line from start to end.
  // The sample points are equally spaced.
  sample_points_.resize(n_sample_points_);
  Double3 direction = end - start;
  direction /= (n_sample_points_ - 1);
  for (size_t i = 0; i < n_sample_points_; i++) {
    sample_points_[i] = start + direction * i;
  }
}

}  // namespace bdm
