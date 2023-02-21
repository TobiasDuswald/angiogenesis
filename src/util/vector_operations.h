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

#ifndef VECTOR_OPERATIONS_H_
#define VECTOR_OPERATIONS_H_

#include <cmath>
#include "biodynamo.h"

namespace bdm {

/// Returns the unit vector defined by phi and theta in spherical coordinates.
/// Helper function for VectorOnConeAroundAxis.
Double3 VectorOnUnitCone(double phi, double theta);

/// Returns a unit vector that lies on the cone defined by (phi, theta) around
/// axis. Intended for behaviors where phi is random.
Double3 VectorOnConeAroundAxis(const Double3& axis, double phi, double theta);

}  // namespace bdm

#endif  // VECTOR_OPERATIONS_H_
