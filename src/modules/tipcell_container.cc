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
#include "core/functor.h"
#include "modules/vessel.h"

namespace bdm {

class AddTipCellsToContainer : public Functor<void, Agent*, AgentHandle> {
 public:
  explicit AddTipCellsToContainer(
      std::vector<std::vector<uint64_t>>* tip_indices, AgentFlatIdxMap* map)
      : tip_indices_(tip_indices), flat_idx_map_(map) {
    ti_ = ThreadInfo::GetInstance();
  }

  void operator()(Agent* agent, AgentHandle ah) {
    auto* vessel = dynamic_cast<Vessel*>(agent);
    if (vessel) {
      if (vessel->GetDaughterLeft() == nullptr) {
        auto idx = flat_idx_map_->GetFlatIdx(ah);
        tip_indices_->operator[](ti_->GetMyThreadId()).push_back(idx);
      }
    }
  }

 private:
  std::vector<std::vector<uint64_t>>* tip_indices_;
  AgentFlatIdxMap* flat_idx_map_;
  ThreadInfo* ti_;
};

class FilterForVessels : public Functor<bool, Agent*> {
 public:
  FilterForVessels() = default;
  bool operator()(Agent* agent) {
    auto* vessel = dynamic_cast<Vessel*>(agent);
    if (vessel) {
      return true;
    } else {
      return false;
    }
  }
};

TipCellContainer::TipCellContainer() {
  rm_ = Simulation::GetActive()->GetResourceManager();
  ti_ = ThreadInfo::GetInstance();
  Update();
}

std::pair<size_t, size_t> TipCellContainer::Get2DIndex(size_t idx) const {
  assert(idx < global_indices_[global_indices_.size() - 1] &&
         "Out of bounds access.");
  // Find the first element in global_indices_ that is larger than idx
  auto it =
      std::upper_bound(global_indices_.begin(), global_indices_.end(), idx);
  // Thread index of the vector that contains the element
  auto thread_id = std::distance(global_indices_.begin(), it);
  // Index of the element in the vector
  auto offset = (thread_id == 0) ? 0 : global_indices_[thread_id - 1];
  auto local_idx = idx - offset;
  return std::make_pair(thread_id, local_idx);
}

void TipCellContainer::Update() {
  num_elements_ = 0;
  tip_indices_.clear();
  tip_indices_.resize(ti_->GetMaxThreads());
  global_indices_.clear();
  global_indices_.resize(ti_->GetMaxThreads());
  flat_idx_map_.Update();
  if (rm_) {
    // Iterate over all agents and add tip cell position to tip_indices_
    AddTipCellsToContainer add_tip_cell_to_container(&tip_indices_,
                                                     &flat_idx_map_);
    FilterForVessels filter_for_vessels;
    rm_->ForEachAgentParallel(add_tip_cell_to_container, &filter_for_vessels);
  }
  // Determine global indices
  global_indices_[0] = tip_indices_[0].size();
  for (size_t i = 1; i < global_indices_.size(); i++) {
    global_indices_[i] = global_indices_[i - 1] + tip_indices_[i].size();
  }
  num_elements_ = global_indices_[global_indices_.size() - 1];

  std::cout << "TipCellContainer::Update() - num_elements_ = " << num_elements_
            << std::endl;
}

size_t TipCellContainer::size() const { return num_elements_; }

const Double3& TipCellContainer::operator[](size_t idx) const {
  auto index = Get2DIndex(idx);

  AgentHandle ah =
      flat_idx_map_.GetAgentHandle(tip_indices_[index.first][index.second]);
  auto* agent = rm_->GetAgent(ah);
  auto* vessel = dynamic_cast<Vessel*>(agent);
  return vessel->GetMassLocation();
}
}  // namespace bdm
