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

#include "biodynamo.h"

namespace bdm {

/// @brief This function registers all collectors for the time series object.
void DefineAndRegisterCollectors();

/// @brief This function visualizes the information collected via the time
/// series.
void PlotAndSaveTimeseries();

}  // namespace bdm

#endif  // ANALYSIS_H_
