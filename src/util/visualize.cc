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

#include <ctime>
#include <iostream>
#include <numeric>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSystem.h"

#include "biodynamo.h"

namespace bdm {

int PlotAndSaveTimeseries(bool save_json) {
  // Get pointers for simulation and TimeSeries data
  auto sim = Simulation::GetActive();
  auto *ts = sim->GetTimeSeries();

  // // Save the TimeSeries Data as JSON to the folder <date_time>
  if (save_json) {
    ts->SaveJson(Concat(sim->GetOutputDir(), "/data.json"));
  }

  // Create a bdm LineGraph that visualizes the TimeSeries data
  bdm::experimental::LineGraph g1(ts, "my result", "Time", "Number of agents",
                                  true);
  g1.Add("quiescent_cells", "quiescent", "L", kBlue, 1.0, 1);
  g1.Add("G1_cells", "G1", "L", kOrange, 1.0, 1);
  g1.Add("SG2_cells", "SG2", "L", kGreen, 1.0, 1);
  g1.Add("hypoxic_cells", "hypoxic", "L", kRed, 1.0, 1);
  g1.Add("dead_cells", "dead", "L", kGray, 1.0, 1);

  g1.Draw();
  g1.SaveAs(Concat(sim->GetOutputDir(), "/cell_timeseries"), {".svg", ".png"});

  // Print info for user to let him/her know where to find simulation results
  std::string info =
      Concat("<PlotAndSaveTimeseries> ", "Results of simulation were saved to ",
             sim->GetOutputDir(), "/");
  std::cout << "Info: " << info << std::endl;

  return 0;
}

}  // namespace bdm
