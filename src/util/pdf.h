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

#ifndef PDF_H
#define PDF_H

namespace bdm {

// PDF for generalized extreme value distribution. The parameter xi is the
// shape parameter. The parameters location and scale transform the distribution
// to F(x, xi) = G((x - location) / scale, xi) / scale.
double gev_pdf(double x, double location, double sigma, double xi);

/// PDF for wald distribution. The parameters location and scale transform the
/// distribution to F(x) = G((x - location) / scale) / scale.
double wald_pdf(double x, double location, double sigma);

}  // namespace bdm

#endif  // PDF_H
