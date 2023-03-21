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

#include <cassert>
#include <cmath>
#include <iostream>

namespace bdm {

double gev_pdf(double x, double location, double scale, double xi) {
  // See
  // https://docs.scipy.org/doc/scipy/reference/generated/
  // scipy.stats.genextreme.html#scipy.stats.genextreme
  double f{0};
  double y = (x - location) / scale;
  if (xi > 0 && y > 1 / xi) {
    // PDF is not defined for y > 1/xi. Return 0 for numerical integration.
    return 0;
  }
  if (xi < 0 && y < 1 / xi) {
    // PDF is not defined for y < 1/xi. Return 0 for numerical integration.
    return 0;
  }
  if (xi == 0) {
    f = std::exp(-std::exp(-y)) * std::exp(-y);
  } else {
    f = std::exp(-std::pow(1 - xi * y, 1 / xi)) *
        std::pow(1 - xi * y, 1 / xi - 1);
  }
  return f / scale;
};

double wald_pdf(double x, double location, double scale) {
  // See
  // https://docs.scipy.org/doc/scipy/reference/generated/
  // scipy.stats.wald.html#scipy.stats.wald
  double f{0};
  double y = (x - location) / scale;
  if (y <= 0) {
    // PDF is not defined for y <= 0. Return 0 for numerical integration.
    return 0;
  }
  f = 1 / std::sqrt(2 * M_PI * y * y * y) *
      std::exp(-(y - 1) * (y - 1) / (2 * y));
  return f / scale;
};

}  // namespace bdm
