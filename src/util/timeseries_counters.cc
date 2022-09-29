// -----------------------------------------------------------------------------
//
// Copyright (C) 2022 CERN, TUM, and UT Austin. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
//
// -----------------------------------------------------------------------------

#include "timeseries_counters.h"
#include "modules/tumor_cell.h"

namespace bdm {

double CountQuiescent(Simulation* sim) {
  auto cond = L2F([](Agent* a) {
    auto* tumor_cell_local = dynamic_cast<TumorCell*>(a);
    if (tumor_cell_local != nullptr) {
      return (tumor_cell_local->GetCellState() == CellState::kQuiescent);
    } else {
      return false;
    }
  });
  return static_cast<double>(bdm::experimental::Count(sim, cond));
};

double CountG1(Simulation* sim) {
  auto cond = L2F([](Agent* a) {
    auto* tumor_cell_local = dynamic_cast<TumorCell*>(a);
    if (tumor_cell_local != nullptr) {
      return (tumor_cell_local->GetCellState() == CellState::kProliferativeG1);
    } else {
      return false;
    }
  });
  return static_cast<double>(bdm::experimental::Count(sim, cond));
};

double CountSG2(Simulation* sim) {
  auto cond = L2F([](Agent* a) {
    auto* tumor_cell_local = dynamic_cast<TumorCell*>(a);
    if (tumor_cell_local != nullptr) {
      return (tumor_cell_local->GetCellState() == CellState::kProliferativeSG2);
    } else {
      return false;
    }
  });
  return static_cast<double>(bdm::experimental::Count(sim, cond));
};

double CountHypoxic(Simulation* sim) {
  auto cond = L2F([](Agent* a) {
    auto* tumor_cell_local = dynamic_cast<TumorCell*>(a);
    if (tumor_cell_local != nullptr) {
      return (tumor_cell_local->GetCellState() == CellState::kHypoxic);
    } else {
      return false;
    }
  });
  return static_cast<double>(bdm::experimental::Count(sim, cond));
};

double CountDead(Simulation* sim) {
  auto cond = L2F([](Agent* a) {
    auto* tumor_cell_local = dynamic_cast<TumorCell*>(a);
    if (tumor_cell_local != nullptr) {
      return (tumor_cell_local->GetCellState() == CellState::kDead);
    } else {
      return false;
    }
  });
  return static_cast<double>(bdm::experimental::Count(sim, cond));
};

double GetSimulatedTime(Simulation* sim) {
  return static_cast<double>(sim->GetScheduler()->GetSimulatedTime());
};

}  // namespace bdm
