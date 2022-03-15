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

#include <gtest/gtest.h>
#include "biodynamo.h"
#include "modules/tumor_cell.h"
#include "util/timeseries_counters.h"

#define TEST_NAME typeid(*this).name()

namespace bdm {

// Creates a Random Cell state that is not equal to state_to_exclude. By
// default, if state_to_exclude=2, the function creates random integers
// 0, 1, 3, ..., CellState::kDead = 4.
int RandomState(CellState state_to_exclude, Random* rng,
                double lower_bound = -1,
                double upper_bound = CellState::kDead) {
  while (true) {
    int random_state = std::ceil(rng->Uniform(lower_bound, upper_bound));
    if (random_state != state_to_exclude) {
      return random_state;
    }
  }
}

// Adds n_state agents in the cell state `state` and n_other agents in other
// cell states to the simulation.
inline void CreateAgents(int n_state, int n_others, CellState state) {
  auto* rm = Simulation::GetActive()->GetResourceManager();
  auto* rnd = Simulation::GetActive()->GetRandom();
  for (int i = 0; i < n_state; i++) {
    auto* tumor_cell = new TumorCell();
    tumor_cell->SetPosition(rnd->UniformArray<3>(10));
    tumor_cell->SetCellState(state);
    rm->AddAgent(tumor_cell);
  }
  for (int i = 0; i < n_state; i++) {
    auto* tumor_cell = new TumorCell();
    tumor_cell->SetPosition(rnd->UniformArray<3>(10));
    tumor_cell->SetCellState(RandomState(state, rnd));
    rm->AddAgent(tumor_cell);
  }
  auto* scheduler = Simulation::GetActive()->GetScheduler();
  scheduler->FinalizeInitialization();
}

// Tests the function CountQuiescent.
TEST(CounterTest, kQuiescent) {
  Simulation simulation(TEST_NAME);
  int n_quiescent = 9;
  int n_others = 7;
  CreateAgents(n_quiescent, n_others, CellState::kQuiescent);
  EXPECT_EQ(n_quiescent, CountQuiescent(&simulation));
}

// Tests the function CountG1.
TEST(CounterTest, kProliferativeG1) {
  Simulation simulation(TEST_NAME);
  int n_proliferativeG1 = 12;
  int n_others = 5;
  CreateAgents(n_proliferativeG1, n_others, CellState::kProliferativeG1);
  EXPECT_EQ(n_proliferativeG1, CountG1(&simulation));
}

// Tests the function CountSG2.
TEST(CounterTest, kProliferativeSG2) {
  Simulation simulation(TEST_NAME);
  int n_proliferativeSG2 = 3;
  int n_others = 18;
  CreateAgents(n_proliferativeSG2, n_others, CellState::kProliferativeSG2);
  EXPECT_EQ(n_proliferativeSG2, CountSG2(&simulation));
}

// Tests the function CountHypoxic.
TEST(CounterTest, kHypoxic) {
  Simulation simulation(TEST_NAME);
  int n_hypoxic = 21;
  int n_others = 4;
  CreateAgents(n_hypoxic, n_others, CellState::kHypoxic);
  EXPECT_EQ(n_hypoxic, CountHypoxic(&simulation));
}

// Tests the function CountDead.
TEST(CounterTest, kDead) {
  Simulation simulation(TEST_NAME);
  int n_dead = 6;
  int n_others = 8;
  CreateAgents(n_dead, n_others, CellState::kDead);
  EXPECT_EQ(n_dead, CountDead(&simulation));
}

}  // namespace bdm
