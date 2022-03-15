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
#include "core/environment/uniform_grid_environment.h"
#include "core/operation/mechanical_forces_op.h"
#include "modules/mechanical_forces.h"
#include "modules/tumor_cell.h"
#include "sim_param.h"

#define TEST_NAME typeid(*this).name()

namespace bdm {

// We implemented a custom force MechanicalInteractionForce and here we test if
// the force produces the correct numeric values.
TEST(ForceTest, AdhesionRepulsion) {
  // Register the simulation parameter
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [](Param* param) {
    auto* sparam = param->Get<SimParam>();
    sparam->viscosity = 2;
    sparam->max_speed = 10.0;
    sparam->adhesion_scale_parameter = 0.0489;
    sparam->repulsive_scale_parameter = 10.0;
  };
  Simulation simulation(TEST_NAME, set_param);
  auto* param = simulation.GetParam();
  auto* sparam = param->Get<SimParam>();
  Double4 computed_force;
  Double3 pos1({0.0, 0.0, 0.0});
  Double3 pos2({0.0, 5.0, 0.0});
  TumorCell t1(pos1, 0);
  TumorCell t2(pos2, 0);
  MechanicalInteractionForce force(sparam->adhesion_scale_parameter,
                                   sparam->repulsive_scale_parameter);
  double expected_y_force{0.0};

  // Set radius, nuclear radius, and action radius.
  t1.SetRadii(1.0, 0.5, 2.0);
  t2.SetRadii(1.0, 0.5, 2.0);

  // Compute force, should be zero since it's far appart right now.
  computed_force = force.Calculate(&t1, &t2);
  EXPECT_EQ(0.0, computed_force[0]);
  EXPECT_EQ(0.0, computed_force[1]);
  EXPECT_EQ(0.0, computed_force[2]);

  // Set radius, nuclear radius, and action radius. We set a large action radius
  // such that we can test if only our adhesion force will be computed.
  t1.SetRadii(0.001, 0.0001, 5.0);
  t2.SetRadii(0.001, 0.0001, 5.0);

  // Compute force and compare to manually calculated force.
  computed_force = force.Calculate(&t1, &t2);
  expected_y_force =
      -0.25 * sparam->adhesion_scale_parameter;  // calculated adhesion force
  EXPECT_EQ(0.0, computed_force[0]);
  EXPECT_EQ((-1) * expected_y_force, computed_force[1]);
  EXPECT_EQ(0.0, computed_force[2]);

  // Set radius, nuclear radius, and action radius. We set a large action radius
  // to zero such that we only consider the repulsion force
  t1.SetRadii(4.0, 2.0, 0.0);
  t2.SetRadii(4.0, 2.0, 0.0);

  // Compute force and compare to manually calculated force.
  computed_force = force.Calculate(&t1, &t2);
  expected_y_force = 9.0 * sparam->repulsive_scale_parameter /
                     64.0;  // calculated repulsion force
  EXPECT_EQ(0.0, computed_force[0]);
  EXPECT_EQ((-1) * expected_y_force, computed_force[1]);
  EXPECT_EQ(0.0, computed_force[2]);

  // Let's move t2 closer to t1 to test the other brach of the if function.
  t2.SetPosition({0.0, 3.0, 0.0});

  // Compute force and compare to manually calculated force.
  computed_force = force.Calculate(&t1, &t2);
  expected_y_force = 7 * sparam->repulsive_scale_parameter /
                     16.0;  // calculated repulsion force
  EXPECT_EQ(0.0, computed_force[0]);
  EXPECT_EQ((-1) * expected_y_force, computed_force[1]);
  EXPECT_EQ(0.0, computed_force[2]);

  // Increase ActionRadius to check comined force
  t1.SetActionRadius(5.0);
  t2.SetActionRadius(5.0);

  // Compute force and compare to manually calculated force.
  computed_force = force.Calculate(&t1, &t2);
  expected_y_force -=
      sparam->adhesion_scale_parameter * 0.49;  // add adhesion force
  EXPECT_EQ(0.0, computed_force[0]);
  EXPECT_EQ((-1) * expected_y_force, computed_force[1]);
  EXPECT_EQ(0.0, computed_force[2]);
}

// The custom force MechanicalInteractionForce is used to compute the
// displacement of the cells. Here we check, if the correct displacements are
// calculated. Note: the numeric force values have been checked above.
TEST(ForceTest, Displacement) {
  // Register the simulation parameter
  Param::RegisterParamGroup(new SimParam());
  auto set_param = [](Param* param) {
    param->simulation_time_step = 0.01;
    param->min_bound = -100;
    param->max_bound = 100;
    auto* sparam = param->Get<SimParam>();
    sparam->viscosity = 2;
    sparam->max_speed = 10.0;
    sparam->adhesion_scale_parameter = 0.0489;
    sparam->repulsive_scale_parameter = 10.0;
    sparam->upper_bound = 100;
    sparam->lower_bound = -100;
  };

  Simulation simulation(TEST_NAME, set_param);
  auto* rm = simulation.GetResourceManager();
  auto* scheduler = simulation.GetScheduler();
  auto* param = simulation.GetParam();
  auto* sparam = param->Get<SimParam>();
  auto* env =
      dynamic_cast<UniformGridEnvironment*>(simulation.GetEnvironment());

  // Set box length manually since diameter_ is not representing the force
  // range. See PR #7.
  env->SetBoxLength(static_cast<int32_t>(std::ceil(2 * 30.5)));

  Double3 computed_displacement;
  Double3 pos1({0.0, 0.0, 0.0});
  Double3 pos2({15.0, 0.0, 0.0});
  Double3 pos3({0.0, 30.0, 0.0});
  Double3 pos4({0.0, 0.0, 7.0});
  TumorCell* t1 = new TumorCell(pos1, 0);
  TumorCell* t2 = new TumorCell(pos2, 0);
  TumorCell* t3 = new TumorCell(pos3, 0);
  TumorCell* t4 = new TumorCell(pos4, 0);
  AgentPointer<TumorCell> agent = t1->GetAgentPtr<TumorCell>();

  // Set radius, nuclear radius, and action radius.
  t1->SetRadii(10, 5, 20);
  t2->SetRadii(10, 5, 20);
  t3->SetRadii(10, 5, 20);
  t4->SetRadii(10, 5, 20);

  // Add agents to simulation
  rm->AddAgent(t1);
  rm->AddAgent(t2);
  rm->AddAgent(t3);
  rm->AddAgent(t4);

  // Use custom force implementation.
  auto* custom_force = new MechanicalInteractionForce(
      sparam->adhesion_scale_parameter, sparam->repulsive_scale_parameter);
  auto* op = scheduler->GetOps("mechanical forces")[0];
  auto* force_implementation = op->GetImplementation<MechanicalForcesOp>();
  force_implementation->SetInteractionForce(custom_force);

  scheduler->FinalizeInitialization();
  scheduler->Simulate(1);

  // Comute displacement
  computed_displacement = agent->GetPosition() - pos1;
  double expected_x{0.003029};
  double expected_y{-1.52812e-05};
  double expected_z{0.02358};
  double tol{0.05};  // accept five percent difference because of col major BDM
  EXPECT_NEAR((-1) * expected_x, computed_displacement[0],
              abs(expected_x * tol));
  EXPECT_NEAR((-1) * expected_y, computed_displacement[1],
              abs(expected_y * tol));
  EXPECT_NEAR((-1) * expected_z, computed_displacement[2],
              abs(expected_z * tol));
}

}  // namespace bdm
