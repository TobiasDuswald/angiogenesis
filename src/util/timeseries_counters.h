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

#ifndef TIMESERIES_COUNTERS_H_
#define TIMESERIES_COUNTERS_H_

#include "biodynamo.h"

namespace bdm {

// Counts the number of TumorCells in CellState::kQuiescent
double CountQuiescent(Simulation* sim);

// Counts the number of TumorCells in CellState::kProliferativeG1
double CountG1(Simulation* sim);

// Counts the number of TumorCells in CellState::kProliferativeSG2
double CountSG2(Simulation* sim);

// Counts the number of TumorCells in CellState::kHypoxic
double CountHypoxic(Simulation* sim);

// Counts the number of TumorCells in CellState::kDead
double CountDead(Simulation* sim);

// Wraps call to bdm::Scheduler::GetSimulatedTime
double GetSimulatedTime(Simulation* sim);

}  // namespace bdm

#endif  // TIMESERIES_COUNTERS_H_
