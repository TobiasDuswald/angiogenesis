// -----------------------------------------------------------------------------
//
// Copyright (C) 2022 CERN, TUM, and UT Austin. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
//
// -----------------------------------------------------------------------------

#ifndef VISUALIZE_H_
#define VISUALIZE_H_

namespace bdm {

// This functions retrieves the collected time series from the active
// simulation, (optionally) saves the results as a JSON file, and plots the
// results. Note that this function assumes certain names given to the different
// timeseries. For details and compatibility, check the .cc implementation.
int PlotAndSaveTimeseries(bool save_json = false);

}  // namespace bdm

#endif  // VISUALIZE_H_
