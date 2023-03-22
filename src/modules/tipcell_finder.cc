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

#include "tipcell_finder.h"
#include <cassert>
#include <chrono>
#include <iostream>
#include "core/util/log.h"
#include "core/util/timing.h"
#include "omp.h"

namespace bdm {

TipCellFinder::TipCellFinder() {
  octree_ = std::make_unique<unibn::Octree<Real3, TipCellContainer>>();
  Update();
}

int TipCellFinder::FindClosestTipCell(const Real3& x) const {
  return octree_->findNeighbor<unibn::L2Distance<Real3>>(x);
}

Real3 TipCellFinder::GetTipCellCenter(int element_id) const {
  return tip_cell_container_[element_id];
}

bool TipCellFinder::IsTipCellInBall(const Real3& x, double r) const {
  if (tip_cell_container_.size() == 0) {
    return false;
  }
  auto idx = FindClosestTipCell(x);
  auto center = GetTipCellCenter(idx);
  auto dist = (center - x).Norm();
  if (dist <= r) {
    return true;
  } else {
    return false;
  }
}

void TipCellFinder::Update() {
  {
    Timing t("TipCellContainer::Update");
    if (update_container_) {
      tip_cell_container_.Update();
    }
  }
  {
    Timing t("Octree::Update");
    if (tip_cell_container_.size() > 0) {
      auto* param = Simulation::GetActive()->GetParam();
      unibn::OctreeParams params;
      params.bucketSize = param->unibn_bucketsize;
      octree_->initialize(tip_cell_container_, params);
    }
  }
  // For all further update calls, we update the container
  update_container_ = true;
}

}  // namespace bdm
