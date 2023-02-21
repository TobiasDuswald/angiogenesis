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

TipCellContainer::TipCellContainer() { Update(); }

void TipCellContainer::Update() {
  auto* rm = Simulation::GetActive()->GetResourceManager();
  assert(rm != nullptr && "ResourceManager is NULL.");
  if (rm) {
    // Iterate over all agents and add tip cell position to
    // element_center_points_
    std::vector<Real3> tmp;
    // ToDo (Tobias): parallelize!
    rm->ForEachAgent([&tmp](Agent* agent) {
      const auto* vessel = dynamic_cast<Vessel*>(agent);
      if (vessel && vessel->GetDaughterLeft() == nullptr) {
        tmp.push_back(vessel->GetMassLocation());
      }
    });
    element_center_points_ = std::move(tmp);
  }
}

size_t TipCellContainer::size() const { return element_center_points_.size(); }

const Double3& TipCellContainer::operator[](size_t idx) const {
  assert(idx < element_center_points_.size() && "Out of bounds access.");
  return element_center_points_[idx];
}
}  // namespace bdm
