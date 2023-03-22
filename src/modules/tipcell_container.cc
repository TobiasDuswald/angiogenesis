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

#include "tipcell_container.h"
#include <cassert>
#include "modules/vessel.h"

namespace bdm {

TipCellContainer::TipCellContainer() {
  rm_ = Simulation::GetActive()->GetResourceManager();
  Update();
}

void TipCellContainer::Update() {
  tip_indices_.clear();
  flat_idx_map_.Update();
  if (rm_) {
    // Iterate over all agents and add tip cell position to
    // element_center_points_
    // std::vector<uint64_t> tmp;
    // ToDo (Tobias): parallelize!
    rm_->ForEachAgent([this](Agent* agent, AgentHandle ah) {
      const auto* vessel = dynamic_cast<Vessel*>(agent);
      if (vessel && vessel->GetDaughterLeft() == nullptr) {
        auto idx = flat_idx_map_.GetFlatIdx(ah);
        tip_indices_.push_back(idx);
      }
    });
    // tip_indices_ = std::move(tmp);
  }
  std::cout << "TipCellContainer::Update() - tip_indices_.size() = "
            << tip_indices_.size() << std::endl;
}

size_t TipCellContainer::size() const { return tip_indices_.size(); }

const Double3& TipCellContainer::operator[](size_t idx) const {
  assert(idx < tip_indices_.size() && "Out of bounds access.");
  AgentHandle ah = flat_idx_map_.GetAgentHandle(tip_indices_[idx]);
  auto* agent = rm_->GetAgent(ah);
  auto* vessel = dynamic_cast<Vessel*>(agent);
  return vessel->GetMassLocation();
}
}  // namespace bdm
