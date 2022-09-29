// --------------------------------------------------------------------------
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
// --------------------------------------------------------------------------

#include "util/analysis.h"
#include "modules/tumor_cell.h"
#include "modules/vessel.h"

namespace bdm {

using experimental::Counter;
using experimental::GenericReducer;

void DefineAndRegisterCollectors() {
  // Get population statistics, i.e. extract data from simulation
  // Get the pointer to the TimeSeries
  auto *ts = Simulation::GetActive()->GetTimeSeries();

  // Define how to get the time values of the TimeSeries
  auto get_time = [](Simulation *sim) {
    return static_cast<double>(sim->GetScheduler()->GetSimulatedTime());
  };

  // Define how to count the agents in quiescent state
  auto is_quiescent = [](Agent *a) {
    auto *cell = dynamic_cast<TumorCell *>(a);
    if (cell) {
      bool tmp = (cell->GetCellState() == CellState::kQuiescent);
      return tmp;
    } else {
      return false;
    }
  };
  ts->AddCollector("q", new Counter<double>(is_quiescent), get_time);

  // Define how to count the agents in state kProliferativeSG2
  auto is_sg2 = [](Agent *a) {
    auto *cell = dynamic_cast<TumorCell *>(a);
    if (cell) {
      bool tmp = (cell->GetCellState() == CellState::kProliferativeSG2);
      return tmp;
    } else {
      return false;
    }
  };
  ts->AddCollector("sg2", new Counter<double>(is_sg2), get_time);

  // Define how to count the agents in state kProliferativeG1
  auto is_g1 = [](Agent *a) {
    auto *cell = dynamic_cast<TumorCell *>(a);
    if (cell) {
      bool tmp = (cell->GetCellState() == CellState::kProliferativeG1);
      return tmp;
    } else {
      return false;
    }
  };
  ts->AddCollector("g1", new Counter<double>(is_g1), get_time);

  // Define how to count the agents in state kHypoxic
  auto is_hypoxic = [](Agent *a) {
    auto *cell = dynamic_cast<TumorCell *>(a);
    if (cell) {
      bool tmp = (cell->GetCellState() == CellState::kHypoxic);
      return tmp;
    } else {
      return false;
    }
  };
  ts->AddCollector("h", new Counter<double>(is_hypoxic), get_time);

  // Define how to count the agents in state kDead
  auto is_dead = [](Agent *a) {
    auto *cell = dynamic_cast<TumorCell *>(a);
    if (cell) {
      bool tmp = (cell->GetCellState() == CellState::kDead);
      return tmp;
    } else {
      return false;
    }
  };
  ts->AddCollector("d", new Counter<double>(is_dead), get_time);

  // Define how to count the agents in state kDead
  auto is_vessel = [](Agent *a) {
    auto *vessel = dynamic_cast<Vessel *>(a);
    if (vessel) {
      return true;
    } else {
      return false;
    }
  };
  ts->AddCollector("v", new Counter<double>(is_vessel), get_time);
}

void PlotAndSaveTimeseries() {
  // Get pointers for simulation and TimeSeries data
  auto sim = Simulation::GetActive();
  auto *ts = sim->GetTimeSeries();

  // Save the TimeSeries Data as JSON to the folder <date_time>
  ts->SaveJson(Concat(sim->GetOutputDir(), "/time-series-data.json"));

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g1(ts, "TumorCell count", "Time",
                                  "Number of agents", true);
  g1.Add("q", "Q", "L", kBlue, 1.0);
  g1.Add("sg2", "SG2", "L", kGreen, 1.0);
  g1.Add("g1", "G1", "L", kOrange, 1.0);
  g1.Add("h", "H", "L", kGray, 1.0);
  g1.Add("d", "D", "L", kBlack, 1.0);
  g1.Draw();
  g1.SaveAs(Concat(sim->GetOutputDir(), "/tumor_cells"), {".svg", ".png"});

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g2(ts, "Vessel count", "Time",
                                  "Number of agents", true);
  g2.Add("v", "Vessel", "L", kBlue, 1.0);
  g2.Draw();
  g2.SaveAs(Concat(sim->GetOutputDir(), "/vessel_agents"), {".svg", ".png"});
}
}  // namespace bdm
