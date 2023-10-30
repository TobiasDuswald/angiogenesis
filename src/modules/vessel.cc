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

#include "vessel.h"
#include "angiogenesis_simulation.h"
#include "modules/tumor_cell.h"
#include "sim_param.h"

namespace bdm {

void Vessel::RunDiscretization() {
  if (!CanGrow() && IsTerminal()) {
    // For vessel agents that are part of the initial vasculature, we do not
    // execute the Discretization() function. The discretization generates new
    // vessels, which are allowed to grow and therefore also secrete
    // nutrients. We do not want our initial vasculature to supply nutrients.
    return;
  }
  // -------------------------------------------------------------------------
  // Modified discretization function from neurite_element.h
  // Rewritten to allow diameter decrease along the vessel
  // -------------------------------------------------------------------------
  if (!IsTerminal()) {
    // if the neurite element is not terminal, we do not split it
    return;
  }
  constexpr double kMaxLength = 10;
  if (GetActualLength() > kMaxLength) {
    auto* new_vessel = SplitVessel(0.1);
    // Set the diameter of the new neurite
    constexpr double kMinDiameter = 5;
    constexpr double kMaxDiameter = 20;
    constexpr double kDiameterDecay = 0.98;
    double diameter = GetDiameter() * kDiameterDecay;
    diameter = std::max(diameter, kMinDiameter);
    diameter = std::min(diameter, kMaxDiameter);
    new_vessel->SetDiameter(diameter);
    new_vessel->UpdateVolume();
  }
}

bool Vessel::IsTipCell() const { return (IsTerminal() && can_grow_); }

bool Vessel::IsStalkCell() const {
  if (IsTerminal()) {
    return false;
  }
  // Vessel must have a left daughter; Check if it is a tip cell
  const auto* daughter = dynamic_cast<const Vessel*>(GetDaughterLeft().Get());
  return daughter->IsTipCell();
}

double Vessel::GetSurfaceArea() const {
  // Vessels are assumed to be cylindrical.
  return Math::kPi * GetDiameter() * GetActualLength();
}

Vessel* Vessel::SplitVessel(real_t distal_portion) {
  neuroscience::SplitNeuriteElementEvent event(distal_portion);
  CreateNewAgents(event, {this});
  // return bdm_static_cast<Vessel*>(event.new_agent[0]);
  return bdm_static_cast<Vessel*>(event.existing_agent);
}

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
  double my_vegf_concentration = dg_guide_->GetValue(dendrite->GetPosition());
  if (my_vegf_concentration < sparam->vegf_threshold_sprouting) {
    return;
  }
  // Get gradient and see if gradient is above a threshold
  Double3 gradient;
  dg_guide_->GetGradient(dendrite->GetPosition(), &gradient);
  if (gradient.Norm() < sparam->vegf_grad_threshold_apical_growth) {
    return;
  }

  /// 2. Check if the next tip cell is at least a given distance away
  auto* asim = static_cast<AngiogenesisSimulation*>(sim);
  auto* finder = asim->GetTipCellFinder();
  if (finder->IsTipCellInBall(dendrite->GetMassLocation(),
                              sparam->min_dist_to_tip_cell)) {
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
  const double sprouting_probability =
      sparam->sprouting_rate * sim->GetParam()->simulation_time_step;
  if (random->Uniform() < sprouting_probability) {
    // Compute sprouting direction on cone around gradient
    double phi = random->Uniform(2 * Math::kPi);
    double theta = random->Uniform(0.25, 0.80);
    Double3 sprouting_direction = VectorOnConeAroundAxis(gradient, phi, theta);

    // Branch vessel
    auto* new_neurite =
        dendrite->Branch(dendrite->GetDiameter(), sprouting_direction,
                         dendrite->GetActualLength() / 2);

    // Set the diameter of the new neurite
    constexpr double kMinDiameter = 5;
    constexpr double kMaxDiameter = 20;
    constexpr double kDiameterDecay = 0.8;

    double new_diameter = dendrite->GetDiameter() * kDiameterDecay;
    new_diameter = std::max(new_diameter, kMinDiameter);
    new_diameter = std::min(new_diameter, kMaxDiameter);
    new_neurite->SetDiameter(new_diameter);
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
  if (gradient.Norm() > 0.016) {
    // This indicates that we're now in the tumor region, stop growth. No longer
    // counted as tip cell.
    auto* vessel = dynamic_cast<Vessel*>(agent);
    vessel->ProhibitGrowth();
  }

  /// 3. If vessel is close to an Tumor cell, we interrupt the growth as well
  double vegf_concentration = dg_guide_->GetValue(dendrite->GetPosition());
  double decision_quotient = std::abs(gradient.Norm() / vegf_concentration);
  if (decision_quotient < sparam->apical_growth_quotient_threshold) {
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

void LineContinuumInteraction::Initialize(const NewAgentEvent& event) {
  Base::Initialize(event);
  const auto* other =
      dynamic_cast<LineContinuumInteraction*>(event.existing_behavior);
  interaction_rate_ = other->interaction_rate_;
}

void LineContinuumInteraction::Run(Agent* agent) {
  auto* vessel = dynamic_cast<Vessel*>(agent);

  // If the behaviour is assigned to a vessel and it is not part of the
  // initial vascular network, we do supply the nutrients.
  if (vessel) {
    // If we secrete and don't consume, then we only consider vessels that
    // can grow.
    if (!vessel->CanGrow()) {
      return;
    }

    // Exclude tip and stalk cells
    if (vessel->IsTipCell() || vessel->IsStalkCell()) {
      return;
    }

    // Get the pointers to the diffusion grids: nutrients, VEGF, DOX, and TRA
    auto* sim = Simulation::GetActive();
    const auto* rm = sim->GetResourceManager();
    auto* param = sim->GetParam();
    const double simulation_time_step = param->simulation_time_step;
    auto* dg_nutrients = bdm_static_cast<DiffusionGrid*>(
        rm->GetContinuum(Substances::kNutrients));
    auto* dg_vegf =
        bdm_static_cast<DiffusionGrid*>(rm->GetContinuum(Substances::kVEGF));
    auto* dg_dox =
        bdm_static_cast<DiffusionGrid*>(rm->GetContinuum(Substances::kDOX));
    auto* dg_tra =
        bdm_static_cast<DiffusionGrid*>(rm->GetContinuum(Substances::kTRA));

    // This is very simple way to supply nutrients. In the current setup,
    // bifurcations add more nutrients than other elements. This should
    // probably be changed. Also note that the behavior makes sense if the
    // vessel is in the same scale as the diffusion grid, or smaller.

    bool skip_first_weight = false;
    const auto* mother_ptr = dynamic_cast<Vessel*>(vessel->GetMother().Get());
    if (mother_ptr &&
        (vessel->GetAgentPtr<Vessel>() == mother_ptr->GetDaughterRight())) {
      skip_first_weight = true;
    }

    // Get start and and point
    const auto& end = vessel->GetMassLocation();
    Double3 start;
    if (mother_ptr) {
      start = mother_ptr->GetMassLocation();
    } else {
      // If mother is soma, start with middle of vessel
      start = vessel->GetPosition();
    }
    double distance = (end - start).Norm();

    // Problem: If the vessel grows, this will not be updated.
    if (!init_) {
      // In the first step, we define how granular the supply is. We do this
      // by demanding that the distance between two sampling points is roughly
      // half the box length of the discretization. Note that we compute the
      // the number of sample points here and not in the constructor because
      // the BoxLength is not initialized in the DiffusionGrid constructor.
      init_ = true;
      std::array<double, 4> box_length = {
          dg_nutrients->GetBoxLength(), dg_vegf->GetBoxLength(),
          dg_dox->GetBoxLength(), dg_tra->GetBoxLength()};
      // Get the minimum box length
      smallest_voxel_size_ =
          *std::min_element(box_length.begin(), box_length.end());
    }
    n_sample_points_ = std::max(
        3,
        static_cast<int>(std::ceil(2 * distance / smallest_voxel_size_ + 1)));
    if (n_sample_points_ != sample_points_.size()) {
      // Update the number of samples and the weights
      ComputeWeights();
    }

    // Compute the sample points
    ComputeSamplePoints(start, end);

    // Modify the continuum values
    std::array<DiffusionGrid*, 4> dg_array = {dg_nutrients, dg_vegf, dg_dox,
                                              dg_tra};
    const double surface = vessel->GetSurfaceArea();
    for (int j = 0; j < 4; j++) {
      auto* dg = dg_array[j];
      const double rate = interaction_rate_[j];
      if (rate == 0) {
        continue;
      }
      for (size_t i = 0; i < n_sample_points_; i++) {
        // We skip the first weight if we are in the right daughter vessel.
        if (skip_first_weight && i == 0) {
          continue;
        }

        // Update the diffusion grid
        double delta_concentration =
            rate * sample_weights_[i] * surface * simulation_time_step;
        dg->ChangeConcentrationBy(sample_points_[i], delta_concentration,
                                  InteractionMode::kLogistic, true);
      }
    }
  }
}

void LineContinuumInteraction::ComputeWeights() {
  sample_weights_.resize(n_sample_points_);
  for (size_t i = 0; i < n_sample_points_; i++) {
    sample_weights_[i] = 1 / static_cast<double>(n_sample_points_ - 1);
  }
  sample_weights_[0] *= 0.5;
  sample_weights_[n_sample_points_ - 1] *= 0.5;
}

void LineContinuumInteraction::ComputeSamplePoints(const Double3& start,
                                                   const Double3& end) {
  // Computes the sample points distributed along the line from start to end.
  // The sample points are equally spaced.
  sample_points_.resize(n_sample_points_);
  Double3 direction = end - start;
  direction /= static_cast<double>(n_sample_points_ - 1);
  for (size_t i = 0; i < n_sample_points_; i++) {
    sample_points_[i] = start + direction * static_cast<double>(i);
  }
}

}  // namespace bdm
