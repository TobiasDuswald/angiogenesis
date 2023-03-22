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

#ifndef TIP_CELL_CONTAINER_H_
#define TIP_CELL_CONTAINER_H_
#include <vector>
#include "biodynamo.h"
#include "core/container/agent_flat_idx_map.h"
#include "modules/vessel.h"

namespace bdm {

class TipCellContainer {
 private:
  std::vector<uint64_t> tip_indices_;
  ResourceManager* rm_;
  AgentFlatIdxMap flat_idx_map_;

 public:
  TipCellContainer();

  /// Copies the element center points of the mesh to element_center_points via
  /// an intermediate vector. Copies in serial.
  void Update();

  /// Returns the number of tip cells.
  size_t size() const;

  /// Retuns the center coordinate of a element in mesh_ labeled by idx.
  const Real3& operator[](size_t idx) const;

  AgentPointer<Vessel> GetAgent(size_t idx) const {
    auto* agent =
        rm_->GetAgent(flat_idx_map_.GetAgentHandle(tip_indices_[idx]));
    Vessel* vessel = dynamic_cast<Vessel*>(agent);
    if (vessel) {
      return vessel->GetAgentPtr<Vessel>();
    } else {
      return AgentPointer<Vessel>();
    }
  }
};

}  // namespace bdm

#endif  // TIP_CELL_CONTAINER_H_
