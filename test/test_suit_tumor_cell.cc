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
// #include "modules/mechanical_forces.h"
#include "sim_param.h"

#define TEST_NAME typeid(*this).name()

namespace bdm {

// This test targets the TumorCell::Divide() member. We test if it actually
// creates a second cell and if the radii are correctly computed for the two
// daughter cells.
TEST(TumorCellTest, CellDivision) {
  // Register the simulation parameter
  Param::RegisterParamGroup(new SimParam());

  auto set_param = [&](Param* param) {
    param->min_bound = -100;
    param->max_bound = 100;
  };

  // Create simulation
  Simulation simulation(TEST_NAME, set_param);
  // Use no force
  auto* scheduler = simulation.GetScheduler();
  scheduler->UnscheduleOp(scheduler->GetOps("mechanical forces")[0]);

  // Add one growing cell to the simulation
  auto* rm = simulation.GetResourceManager();
  auto* tumor_cell = new TumorCell({0, 0, 0}, 0);
  double radius = 5.0;
  double volume = 4.0 * Math::kPi * pow(radius, 3) / 3.0;
  double displacement_scale_factor{4.0};
  double eps = 1e-3;
  double volume_ratio_max{1.1};
  tumor_cell->SetVolume(volume);
  tumor_cell->SetRadius(radius);
  tumor_cell->SetDisplacementScaleFactor(displacement_scale_factor);
  rm->AddAgent(tumor_cell);
  AgentPointer<TumorCell> tumor_cell_ptr = tumor_cell->GetAgentPtr<TumorCell>();

  // Make cell divide
  tumor_cell_ptr->Divide();
  scheduler->Simulate(1);

  // Test if we have two cells.
  EXPECT_EQ(2u, rm->GetNumAgents());

  // Test if volumes of the split up cells add up.
  double total_volume{0.0};
  rm->ForEachAgent([&](Agent* agent) {
    TumorCell* tc = dynamic_cast<TumorCell*>(agent);
    total_volume += tc->GetVolume();
  });
  EXPECT_FLOAT_EQ(volume, total_volume);

  // Test if cells are no longer located at origin but also not too far off.
  // More precise, if their center is on the right shell on a sphere.
  rm->ForEachAgent([&](Agent* agent) {
    TumorCell* tc = dynamic_cast<TumorCell*>(agent);
    double displacement = tc->GetPosition().Norm();
    EXPECT_LT(radius / displacement_scale_factor / (1 + volume_ratio_max) + eps,
              displacement);
    EXPECT_GT(radius / displacement_scale_factor / (2 + volume_ratio_max) *
                      (1 + volume_ratio_max) -
                  eps,
              displacement);
  });

  // Test if radius and diameter are correctly coupled
  rm->ForEachAgent([&](Agent* agent) {
    TumorCell* tc = dynamic_cast<TumorCell*>(agent);
    EXPECT_FLOAT_EQ(tc->GetRadius() * 2, tc->GetDiameter());
    EXPECT_NE(tc->GetActionRadius() * 2, tc->GetDiameter());
  });
}

// This test targets TumorCell::ChangeVolume(double) to see if it changes
// the volume and the radius correctly.
TEST(TumorCellTest, ChangeVolume) {
  auto set_param = [&](Param* param) { param->simulation_time_step = 0.01; };
  // Create simulation
  Simulation simulation(TEST_NAME, set_param);

  // Add one growing cell to the simulation
  auto* rm = simulation.GetResourceManager();
  auto* tumor_cell = new TumorCell({0, 0, 0}, 0);
  AgentPointer<TumorCell> tumor_cell_ptr = tumor_cell->GetAgentPtr<TumorCell>();
  rm->AddAgent(tumor_cell);

  //
  double radius = 5.0;
  double volume = 4.0 * Math::kPi * pow(radius, 3) / 3.0;
  double new_radius = 6.0;
  double new_volume = 4.0 * Math::kPi * pow(new_radius, 3) / 3.0;
  tumor_cell_ptr->SetVolume(volume);
  tumor_cell_ptr->SetRadius(radius);

  // Get cell volume at the beginning
  EXPECT_EQ(volume, tumor_cell_ptr->GetVolume());
  // 100 to compenstate ts 0.01
  tumor_cell_ptr->ChangeVolume((new_volume - volume) * 100);

  // Get cell volume at the beginning after 3 timesteps
  EXPECT_DOUBLE_EQ(new_volume, tumor_cell_ptr->GetVolume());
  EXPECT_DOUBLE_EQ(new_radius, tumor_cell_ptr->GetRadius());
}

// Test if Apoptosis decreases the cell volume correctly. The decrease depends
// on the targeted size, the current size, and the apoptosis_duration.
TEST(TumorCellTest, ApoptosisVolumeDecrease) {
  auto set_param = [&](Param* param) { param->simulation_time_step = 0.01; };
  // Create simulation
  Simulation simulation(TEST_NAME, set_param);

  // Add one growing cell to the simulation
  auto* rm = simulation.GetResourceManager();
  auto* tumor_cell = new TumorCell({0, 0, 0}, 0);
  AgentPointer<TumorCell> tumor_cell_ptr = tumor_cell->GetAgentPtr<TumorCell>();
  rm->AddAgent(tumor_cell);

  // Compute the apoptosis duration for the cell.
  double radius = 5.0;
  double nuclear_radius = 2.0;
  double apoptosis_duration = 1000;
  tumor_cell_ptr->SetRadius(radius);
  tumor_cell_ptr->SetNuclearRadius(nuclear_radius);
  tumor_cell_ptr->ComputeApoptosisVolumeDecrease(apoptosis_duration);

  // Change Volume multiple times as in simulation.
  for (int i = 0;
       i < static_cast<int>(apoptosis_duration /
                            simulation.GetParam()->simulation_time_step) +
               1;
       i++) {
    tumor_cell_ptr->ChangeVolume(tumor_cell_ptr->GetGrowthRate());
  }

  // Expect volume to be smaller than the nuclear volume. In the simulation, we
  // then remove the agent.
  EXPECT_GT(4.0 / 3.0 * Math::kPi * pow(nuclear_radius, 3),
            tumor_cell_ptr->GetVolume());
}

// Test UpdateHypoxic behavior.
TEST(TumorCellTest, HypoxicTransition) {
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [&](Param* param) {
    param->simulation_time_step = 0.01;
    param->Get<SimParam>()->hypoxic_threshold = 0.3;
  };
  // Create simulation
  Simulation simulation(TEST_NAME, set_param);

  // Define substance
  int substance_id = 0;
  std::string substance_name = "VEGF";
  ModelInitializer::DefineSubstance(substance_id, substance_name, 0, 0, 20);
  auto SetInitialValuesGridVEGF = [&](double x, double y, double z) {
    if (z > 5) {
      return 0.5;
    } else {
      return 0.0;
    }
  };
  ModelInitializer::InitializeSubstance(substance_id, SetInitialValuesGridVEGF);

  // Add cells
  Double3 pos_hypoxic = {0, 0, 0};
  Double3 pos_quiescent = {0, 0, 10};
  auto* rm = simulation.GetResourceManager();
  auto* tumor_cell_1 = new TumorCell(pos_hypoxic, CellState::kQuiescent);
  tumor_cell_1->AddBehavior(new UpdateHypoxic(substance_id));
  auto* tumor_cell_2 = new TumorCell(pos_quiescent, CellState::kHypoxic);
  tumor_cell_2->AddBehavior(new UpdateHypoxic(substance_id));
  rm->AddAgent(tumor_cell_1);
  rm->AddAgent(tumor_cell_2);

  // Simulate 1 iterations
  auto* scheduler = simulation.GetScheduler();
  scheduler->UnscheduleOp(scheduler->GetOps("mechanical forces")[0]);
  scheduler->FinalizeInitialization();
  scheduler->Simulate(1);

  EXPECT_EQ(CellState::kHypoxic, tumor_cell_1->GetCellState());
  EXPECT_EQ(CellState::kQuiescent, tumor_cell_2->GetCellState());
}

// Test if Cells only secrete if they are hypoxic
TEST(TumorCellTest, HypoxicSecretion) {
  auto set_param = [&](Param* param) { param->simulation_time_step = 0.01; };
  // Create simulation
  Simulation simulation(TEST_NAME, set_param);

  // Define substance
  std::string substance_name = "VEGF";
  ModelInitializer::DefineSubstance(0, substance_name, 0, 0, 20);
  auto SetInitialValuesGridVEGF = [&](double x, double y, double z) {
    return 0.0;
  };
  ModelInitializer::InitializeSubstance(0, SetInitialValuesGridVEGF);

  // Add cells
  Double3 pos_hypoxic = {0, 0, 0};
  Double3 pos_quiescent = {0, 0, 10};
  auto* rm = simulation.GetResourceManager();
  auto* tumor_cell_1 = new TumorCell(pos_hypoxic, CellState::kHypoxic);
  tumor_cell_1->AddBehavior(new HypoxicSecretion(substance_name, 1.0));
  auto* tumor_cell_2 = new TumorCell(pos_quiescent, CellState::kQuiescent);
  tumor_cell_2->AddBehavior(new HypoxicSecretion(substance_name, 1.0));
  rm->AddAgent(tumor_cell_1);
  rm->AddAgent(tumor_cell_2);

  // Simulate 10 iterations
  auto* scheduler = simulation.GetScheduler();
  scheduler->UnscheduleOp(scheduler->GetOps("mechanical forces")[0]);
  scheduler->FinalizeInitialization();
  scheduler->Simulate(10);

  auto* dgrid = Simulation::GetActive()->GetResourceManager()->GetDiffusionGrid(
      substance_name);
  auto value_hypoxic = dgrid->GetValue(pos_hypoxic);
  auto value_quiescent = dgrid->GetValue(pos_quiescent);
  EXPECT_EQ(10, value_hypoxic);
  EXPECT_EQ(0, value_quiescent);
}

}  // namespace bdm
