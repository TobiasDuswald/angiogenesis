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

#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include <string>
#include <vector>
#include "biodynamo.h"

namespace bdm {

// Operation to track the min, max and average of all continuum concentrations.
class VerifyContinuum : public StandaloneOperationImpl {
  BDM_OP_HEADER(VerifyContinuum);

  // Apply operation
  void operator()() override;

  // Return the results, e.g. min, max, avg of all concentrations
  std::map<std::string, std::vector<real_t>>* GetResults() {
    return &results_;
  };

  // Get the time points at which the results were recorded
  std::vector<real_t>* GetSimulatedTime() { return &simulated_time_; };

 private:
  // Initialize the results map with the correct keys
  void Initialize();

  // Time points at which the results were recorded
  std::vector<real_t> simulated_time_;

  // Results map
  std::map<std::string, std::vector<real_t>> results_;

  // Flag to indicate if the results map has been initialized
  bool initialized_ = false;
};

/// @brief This function registers all collectors for the time series object.
void DefineAndRegisterCollectors();

/// @brief This function visualizes the information collected via the time
/// series.
void PlotAndSaveTimeseries();

}  // namespace bdm

#endif  // ANALYSIS_H_
