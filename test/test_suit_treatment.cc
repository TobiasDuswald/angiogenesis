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
#include <array>
#include <fstream>
#include "modules/treatment.h"

namespace bdm {

double DaysToMinutes(double days) { return days * 24 * 60; }

TEST(Treatment, Schedule) {
  Treatment treatment;

  // Time in days
  double time = 0;

  // Treatment parameters. Start end in minutes.
  double tra_start_1_ = 102 * 60 * 24;
  double tra_end_1_ = 103 * 60 * 24;
  double tra_start_2_ = 105 * 60 * 24;
  double tra_end_2_ = 106 * 60 * 24;
  double dox_start_ = 108 * 60 * 24;
  double dox_end_ = 109 * 60 * 24;

  // Iterate form 0 to 200 days in steps of 0.1 days and check if the treatment
  // is applied.
  for (int i = 0; i < 2000; i++) {
    time = i * 0.1 * 24 * 60;
    if (time >= tra_start_1_ && time < tra_end_1_) {
      EXPECT_TRUE(treatment.IsTraApplied(time));
    } else if (time >= tra_start_2_ && time < tra_end_2_) {
      EXPECT_TRUE(treatment.IsTraApplied(time));
    } else {
      EXPECT_FALSE(treatment.IsTraApplied(time));
    }

    if (time >= dox_start_ && time < dox_end_) {
      EXPECT_TRUE(treatment.IsDoxApplied(time));
    } else {
      EXPECT_FALSE(treatment.IsDoxApplied(time));
    }
  }
}

TEST(Treatment, VesselPermeabilityODE) {
  Treatment treatment;

  // Time in days
  double time = 0;

  // Treatment parameters. Start end in minutes.
  double tra_start_1_ = 102 * 60 * 24;
  double tra_end_1_ = 103 * 60 * 24;
  double tra_start_2_ = 105 * 60 * 24;
  double tra_end_2_ = 106 * 60 * 24;

  // Iterate form 0 to 200 days in steps of 0.1 days and check if the treatment
  // is applied.
  for (int i = 0; i < 2000; i++) {
    time = i * 0.1 * 24 * 60;
    // Check if The VesselPermeabilityODE is larger than zero if the TRA
    // treatment is applied, smaller than zero else.
    if (time >= tra_start_1_ && time < tra_end_1_) {
      EXPECT_GT(treatment.VesselPermeabilityODE(0.5, time), 0);
    } else if (time >= tra_start_2_ && time < tra_end_2_) {
      EXPECT_GT(treatment.VesselPermeabilityODE(0.5, time), 0);
    } else {
      EXPECT_LT(treatment.VesselPermeabilityODE(0.5, time), 0);
    }
  }
}

TEST(Treatment, PrecomputeVesselPermeability) {
  Treatment treatment;

  // Treatment parameters. Start end in minutes.
  double tra_start_1_ = 102 * 60 * 24;
  double tra_end_1_ = 103 * 60 * 24;
  double tra_start_2_ = 105 * 60 * 24;
  double tra_end_2_ = 106 * 60 * 24;

  // Time in days to precompute the vessel permeability
  double t_end = 200;
  // Time step in days
  double dt = 0.1;
  // Time step for the ODE solver
  double dt_ode = 0.0001;

  // Precompute the vessel permeability
  treatment.PrecomputeVesselPermeability(
      DaysToMinutes(t_end), DaysToMinutes(dt), DaysToMinutes(dt_ode));
  auto& vessel_permeability = treatment.GetVesselPermeability();

  // Iterate over the precomputed vessel permeability and check if the
  // permeability increases if the TRA treatment is applied, and decreases
  // otherwise.
  for (int i = 1021; i < vessel_permeability.size(); i++) {
    double t = i * dt * 24 * 60;
    std::cout << i << " " << t << " " << vessel_permeability[i] << std::endl;
    if (t > tra_start_1_ && t <= tra_end_1_) {
      EXPECT_GT(vessel_permeability[i] - vessel_permeability[i - 1], 0);
    } else if (t > tra_start_2_ && t <= tra_end_2_) {
      EXPECT_GT(vessel_permeability[i] - vessel_permeability[i - 1], 0);
    } else {
      EXPECT_LT(vessel_permeability[i] - vessel_permeability[i - 1], 0);
    }
  }

  // Write the vector to a txt file
  std::ofstream file;
  file.open("vessel_permeability.txt");
  for (int i = 0; i < vessel_permeability.size(); i++) {
    file << vessel_permeability[i] << std::endl;
  }
}

}  // namespace bdm
