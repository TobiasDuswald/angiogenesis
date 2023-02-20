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
#include "sim_param.h"

namespace bdm {

using experimental::Counter;
using experimental::GenericReducer;

void VerifyContinuum::operator()() {
  if (!initialized_) {
    // Initialize the results map with the correct keys
    Initialize();
  }

  // Get the current simulation time
  auto *sim = Simulation::GetActive();
  auto *rm = sim->GetResourceManager();
  simulated_time_.push_back(sim->GetScheduler()->GetSimulatedTime());

  // Iterate over all continuum models of type diffusion grid
  rm->ForEachDiffusionGrid([&](DiffusionGrid *grid) {
    // Compute min, max and avg of the current grid
    auto &cn = grid->GetContinuumName();
    const real_t *data = grid->GetAllConcentrations();
    int num_elements = grid->GetNumBoxes();
    real_t min = data[0];
    real_t max = data[0];
    real_t sum = 0;
    for (int i = 0; i < num_elements; i++) {
      min = std::min(min, data[i]);
      max = std::max(max, data[i]);
      sum += data[i];
    }
    // Add the results to the results map
    results_[cn + "_min"].push_back(min);
    results_[cn + "_max"].push_back(max);
    results_[cn + "_avg"].push_back(sum / num_elements);
  });
};

void VerifyContinuum::Initialize() {
  initialized_ = true;
  auto *sim = Simulation::GetActive();
  auto *rm = sim->GetResourceManager();
  // Initialize a map with the correct keys, e.g "<continuum>_min",
  // "<continuum>_max", "<continuum>_avg"
  rm->ForEachDiffusionGrid([&](DiffusionGrid *grid) {
    auto &cn = grid->GetContinuumName();
    results_[cn + "_min"] = std::vector<real_t>();
    results_[cn + "_max"] = std::vector<real_t>();
    results_[cn + "_avg"] = std::vector<real_t>();
  });
};

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
  auto *sparam = sim->GetParam()->Get<SimParam>();
  auto *ts = sim->GetTimeSeries();
  auto *scheduler = sim->GetScheduler();

  // Save the TimeSeries Data as JSON to the folder <date_time>
  ts->SaveJson(Concat(sim->GetOutputDir(), "/time-series-data.json"));

  // Create a bdm LineGraph that visualizes the TimeSeries data
  {
    bdm::experimental::LineGraph g1(ts, "TumorCell count", "Time [min]",
                                    "Number of agents", true);
    g1.Add("q", "Q", "L", kOrange, 1.0);
    g1.Add("sg2", "SG2", "L", kGreen + 2, 1.0);
    g1.Add("g1", "G1", "L", kGreen, 1.0);
    g1.Add("h", "H", "L", kGray + 1, 1.0);
    g1.Add("d", "D", "L", kBlack, 1.0);
    g1.Draw();
    g1.SaveAs(Concat(sim->GetOutputDir(), "/tumor_cells"), {".pdf", ".png"});
  }

  // Transform the data of the timeseries by scaling the x-axis from minutes to
  // days
  constexpr double min_to_days = 1.0 / 60.0 / 24.0;
  bdm::experimental::LinearTransformer lt;
  lt.SetXSlope(min_to_days);
  ts->AddTransformedData("q", "q_days", lt);
  ts->AddTransformedData("sg2", "sg2_days", lt);
  ts->AddTransformedData("g1", "g1_days", lt);
  ts->AddTransformedData("h", "h_days", lt);
  ts->AddTransformedData("d", "d_days", lt);

  // Create line graph that visualizes the TimeSeries data in days
  {
    bdm::experimental::LineGraph g1(ts, "TumorCell count", "Time [days]",
                                    "Number of agents", true);
    g1.Add("q_days", "Q", "L", kOrange, 1.0);
    g1.Add("sg2_days", "SG2", "L", kGreen + 2, 1.0);
    g1.Add("g1_days", "G1", "L", kGreen, 1.0);
    g1.Add("h_days", "H", "L", kGray + 1, 1.0);
    g1.Add("d_days", "D", "L", kBlack, 1.0);
    g1.Draw();
    g1.SaveAs(Concat(sim->GetOutputDir(), "/tumor_cells_days"),
              {".pdf", ".png"});
  }

  // Tranform the data of the timeseries to volume by multiplying the number of
  // agents with the volume of a single cell
  constexpr double volume = 4.0 / 3.0 * M_PI * 0.01 * 0.01 * 0.01;  // mm^3
  bdm::experimental::LinearTransformer lt2;
  lt2.SetXSlope(min_to_days);
  lt2.SetYSlope(volume);
  ts->AddTransformedData("q", "q_volume", lt2);
  ts->AddTransformedData("sg2", "sg2_volume", lt2);
  ts->AddTransformedData("g1", "g1_volume", lt2);
  ts->AddTransformedData("h", "h_volume", lt2);
  ts->AddTransformedData("d", "d_volume", lt2);

  // Create line graph that visualizes the TimeSeries data in volume
  {
    bdm::experimental::LineGraph g1(ts, "Tumor volume", "Time [days]",
                                    "Volume [mm^3]", true);
    g1.Add("q_volume", "Q", "L", kOrange, 1.0);
    g1.Add("sg2_volume", "SG2", "L", kGreen + 2, 1.0);
    g1.Add("g1_volume", "G1", "L", kGreen, 1.0);
    g1.Add("h_volume", "H", "L", kGray + 1, 1.0);
    g1.Add("d_volume", "D", "L", kBlack, 1.0);
    g1.Draw();
    g1.SaveAs(Concat(sim->GetOutputDir(), "/tumor_cells_volume"),
              {".pdf", ".png"});
  }

  // Create a bdm LineGraph that visualizes the TimeSeries data
  {
    bdm::experimental::LineGraph g2(ts, "Vessel count", "Time",
                                    "Number of agents", true);
    g2.Add("v", "Vessel", "L", kBlue, 1.0);
    g2.Draw();
    g2.SaveAs(Concat(sim->GetOutputDir(), "/vessel_agents"), {".pdf", ".png"});
  }

  // Add the TimeSeries from the continuum verification to the TimeSeries
  if (sparam->verify_continuum_values) {
    auto *op = scheduler->GetOps("VerifyContinuum")[0];
    auto *results_continuum =
        op->GetImplementation<VerifyContinuum>()->GetResults();
    auto *sim_time =
        op->GetImplementation<VerifyContinuum>()->GetSimulatedTime();
    // Add all key values pairs to the TimeSeries
    for (auto &pair : *results_continuum) {
      ts->Add(pair.first, *sim_time, pair.second);
    }

    // Create a bdm LineGraph that visualizes the TimeSeries data for the
    // continuum
    {
      // Plot for nutrients
      bdm::experimental::LineGraph g3(ts, "Nutrients", "Time", "Value", true);
      g3.Add("Nutrients_avg", "Nutrients (avg)", "L", kRed, 1.0);
      g3.Add("Nutrients_min", "Nutrients (min)", "L", kGreen, 1.0);
      g3.Add("Nutrients_max", "Nutrients (max)", "L", kBlack, 1.0);
      g3.Draw();
      g3.SaveAs(Concat(sim->GetOutputDir(), "/continuum_values_nutrients"),
                {".pdf", ".png"});
    }
    {
      // Plot for VEGF
      bdm::experimental::LineGraph g3(ts, "VEGF", "Time", "Value", true);
      g3.Add("VEGF_avg", "VEGF (avg)", "L", kRed, 1.0);
      g3.Add("VEGF_min", "VEGF (min)", "L", kGreen, 1.0);
      g3.Add("VEGF_max", "VEGF (max)", "L", kBlack, 1.0);
      g3.Draw();
      g3.SaveAs(Concat(sim->GetOutputDir(), "/continuum_values_vegf"),
                {".pdf", ".png"});
    }
    {
      // Plot for DOX
      bdm::experimental::LineGraph g3(ts, "DOX", "Time", "Value", true);
      g3.Add("DOX_avg", "DOX (avg)", "L", kRed, 1.0);
      g3.Add("DOX_min", "DOX (min)", "L", kGreen, 1.0);
      g3.Add("DOX_max", "DOX (max)", "L", kBlack, 1.0);
      g3.Draw();
      g3.SaveAs(Concat(sim->GetOutputDir(), "/continuum_values_dox"),
                {".pdf", ".png"});
    }
    {
      // Plot for TRA
      bdm::experimental::LineGraph g3(ts, "TRA", "Time", "Value", true);
      g3.Add("TRA_avg", "TRA (avg)", "L", kRed, 1.0);
      g3.Add("TRA_min", "TRA (min)", "L", kGreen, 1.0);
      g3.Add("TRA_max", "TRA (max)", "L", kBlack, 1.0);
      g3.Draw();
      g3.SaveAs(Concat(sim->GetOutputDir(), "/continuum_values_tra"),
                {".pdf", ".png"});
    }
  }
}
}  // namespace bdm
