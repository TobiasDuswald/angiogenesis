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

#ifndef ANGIOGENESIS_SIMULATION_H_
#define ANGIOGENESIS_SIMULATION_H_

#include "biodynamo.h"

namespace bdm {

// Function to create a TumorCell with random properties. This function is
// passed on to the ModelInitializer creating random positions according to some
// distribution.
auto CreateTumorCell(const Double3& position);

// Wrapper to multiple call to CreateTumorCell.
void PlaceTumorCells();

// This function contains the core simulation code. It creates the agents in
// an environment and simulates the system for multiple timesteps.
int Simulate(int argc, const char** argv);

}  // namespace bdm

#endif  // ANGIOGENESIS_SIMULATION_H_
