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
#include "modules/element_container.h"
#include "modules/element_finder.h"
#include "modules/vessel.h"

#define TEST_NAME typeid(*this).name()

namespace bdm {

// --------------------------------------------------------------------------
// Helper function taken from angiogenesis_simulation.cc to construct the
// simulation environment, e.g. the vessels.
// --------------------------------------------------------------------------
void inline PlaceVessel(Double3 start, Double3 end, double compartment_length) {
  auto* rm = Simulation::GetActive()->GetResourceManager();

  // Compute parameters for straight line between start and end.
  Double3 direction = end - start;
  double distance = direction.Norm();
  direction.Normalize();
  int n_compartments = std::floor(distance / compartment_length);

  // Warn if chosen parameters are not selected ideally
  if (abs(n_compartments * compartment_length - distance) > 1e-2) {
    Log::Warning("PlaceVessel", "Vessel will be shorter than expected.");
  }

  // The setup requires us to define a NeuronSoma, which is kind of a left over
  // from the neuroscience module.
  const Double3 tmp = start;
  auto* soma = new neuroscience::NeuronSoma(tmp);
  rm->AddAgent(soma);

  // Define a first neurite
  Vessel v;  // Used for prototype argument (virtual+template not supported c++)
  auto* vessel_compartment_1 =
      dynamic_cast<Vessel*>(soma->ExtendNewNeurite(direction, &v));
  vessel_compartment_1->SetPosition(start +
                                    direction * compartment_length * 0.5);
  vessel_compartment_1->SetMassLocation(start + direction * compartment_length);
  vessel_compartment_1->SetActualLength(compartment_length);
  vessel_compartment_1->SetDiameter(15);
  vessel_compartment_1->ProhibitGrowth();

  Vessel* vessel_compartment_2{nullptr};
  for (int i = 1; i < n_compartments; i++) {
    // Compute location of next vessel element
    Double3 agent_position =
        start + direction * compartment_length * (static_cast<double>(i) + 0.5);
    Double3 agent_end_position =
        start + direction * compartment_length * (static_cast<double>(i) + 1.0);
    // Create new vessel
    vessel_compartment_2 = new Vessel();
    // Set position an length
    vessel_compartment_2->SetPosition(agent_position);
    vessel_compartment_2->SetMassLocation(agent_end_position);
    vessel_compartment_2->SetActualLength(compartment_length);
    vessel_compartment_2->SetRestingLength(compartment_length);
    vessel_compartment_2->SetSpringAxis(direction);
    vessel_compartment_2->SetDiameter(15);
    vessel_compartment_2->ProhibitGrowth();
    // Add Agent to the resource manager
    rm->AddAgent(vessel_compartment_2);
    // Connect vessels (AgentPtr API is currently bounded to base
    // classes but this is a 'cosmetic' problem)
    vessel_compartment_1->SetDaughterLeft(
        vessel_compartment_2->GetAgentPtr<neuroscience::NeuriteElement>());
    vessel_compartment_2->SetMother(
        vessel_compartment_1->GetAgentPtr<neuroscience::NeuronOrNeurite>());
    std::swap(vessel_compartment_1, vessel_compartment_2);
  }
}

TEST(TipCellFinder, Container) {
  neuroscience::InitModule();
  Simulation simulation(TEST_NAME);
  PlaceVessel({0, 0, 0}, {100, 0, 0}, 1);
  PlaceVessel({0, 100, 0}, {100, 100, 0}, 1);
  PlaceVessel({0, 100, 100}, {100, 100, 100}, 1);
  auto* scheduler = Simulation::GetActive()->GetScheduler();
  scheduler->Simulate(1);

  TipCellContainer container;
  // Check the number of tip cells
  EXPECT_EQ(container.size(), 3u);
  // Check the location of the tip cells
  for (size_t i = 0; i < container.size(); i++) {
    EXPECT_LT(container[i][0] - 104, 0);
    EXPECT_GT(container[i][0] - 102, 0);
  }
}

TEST(TipCellFinder, Finder) {
  neuroscience::InitModule();
  Simulation simulation(TEST_NAME);
  const Real3 tip1 = {103, 0, 0};
  const Real3 tip2 = {103, 100, 0};
  PlaceVessel({0, 0, 0}, {100, 0, 0}, 1);
  PlaceVessel({0, 100, 0}, {100, 100, 0}, 1);
  auto* scheduler = Simulation::GetActive()->GetScheduler();
  scheduler->Simulate(1);

  TipCellFinder finder;
  // Check the number of tip cells
  EXPECT_EQ(finder.GetNumberOfTipCells(), 2u);
  // Check the location of the tip cells
  for (size_t i = 0; i < finder.GetNumberOfTipCells(); i++) {
    auto tip = finder.GetTipCellCenter(i);
    bool found_tip = ((tip - tip1).Norm() < 1e-6 || (tip - tip2).Norm() < 1e-6);
    EXPECT_TRUE(found_tip);
  }
  // Check the IsTipCellInBall function
  const Real3 test_point_1 = {103, 0, 0};
  const Real3 test_point_2 = {103, 100, 0};
  const Real3 test_point_3 = {103, 50, 0};
  const Real3 test_point_4 = {103, 30, 0};
  const Real3 test_point_5 = {103, 70, 0};
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_1, 1));
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_2, 1));
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_3, 50 + 1e-6));
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_4, 30 + 1e-6));
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_5, 30 + 1e-6));
  EXPECT_FALSE(finder.IsTipCellInBall(test_point_3, 49));
  EXPECT_FALSE(finder.IsTipCellInBall(test_point_4, 29));
  EXPECT_FALSE(finder.IsTipCellInBall(test_point_5, 29));

  // Add vessel and update tip cell finder
  const Real3 tip3 = {106, 0, 0};
  const Real3 tip4 = {103, 50, 0};  // new vessel hasn't moved
  const Real3 tip5 = {106, 100, 0};
  PlaceVessel({0, 50, 0}, {100, 50, 0}, 1);
  scheduler->Simulate(1);
  finder.Update();

  // Check the number of tip cells
  EXPECT_EQ(finder.GetNumberOfTipCells(), 3u);
  // Check the location of the tip cells
  for (size_t i = 0; i < finder.GetNumberOfTipCells(); i++) {
    auto tip = finder.GetTipCellCenter(i);
    std::cout << tip << std::endl;
    bool found_tip = ((tip - tip3).Norm() < 1e-6 ||
                      (tip - tip4).Norm() < 1e-6 || (tip - tip5).Norm() < 1e-6);
    EXPECT_TRUE(found_tip);
  }
  const Real3 test_point_6 = {106, 0, 0};
  const Real3 test_point_7 = {106, 100, 0};
  const Real3 test_point_8 = {103, 50, 0};
  const Real3 test_point_9 = {103, 30, 0};
  const Real3 test_point_10 = {103, 70, 0};
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_6, 1));
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_7, 1));
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_8, 1));
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_9, 20 + 1e-6));
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_10, 20 + 1e-6));
  EXPECT_TRUE(finder.IsTipCellInBall(test_point_8, 0));  // Self find
  EXPECT_FALSE(finder.IsTipCellInBall(test_point_9, 19));
  EXPECT_FALSE(finder.IsTipCellInBall(test_point_10, 19));
}

}  // namespace bdm
