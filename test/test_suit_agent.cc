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
TEST(AgentTest, CellDivision) {
  // Register the simulation parameter
  Param::RegisterParamGroup(new SimParam());

  auto set_param = [&](Param* param) {
    param->min_bound = -100;
    param->max_bound = 100;
    auto* sparam = param->Get<SimParam>();
    sparam->upper_bound = 100;
    sparam->lower_bound = -100;
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
    TumorCell* tc = bdm_static_cast<TumorCell*>(agent);
    total_volume += tc->GetVolume();
  });
  EXPECT_FLOAT_EQ(volume, total_volume);

  // Test if cells are no longer located at origin but also not too far off.
  // More precise, if their center is on the right shell on a sphere.
  rm->ForEachAgent([&](Agent* agent) {
    TumorCell* tc = bdm_static_cast<TumorCell*>(agent);
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
    TumorCell* tc = bdm_static_cast<TumorCell*>(agent);
    EXPECT_FLOAT_EQ(tc->GetRadius() * 2, tc->GetDiameter());
    EXPECT_NE(tc->GetActionRadius() * 2, tc->GetDiameter());
  });
}

// This test targets TumorCell::ChangeVolume(double) to see if it changes
// the volume and the radius correctly.
TEST(AgentTest, ChangeVolume) {
  // Create simulation
  Simulation simulation(TEST_NAME);

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
TEST(AgentTest, ApoptosisVolumeDecrease) {
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

}  // namespace bdm
