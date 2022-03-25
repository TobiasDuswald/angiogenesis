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

#include "vector_operations.h"
#include <cmath>

namespace bdm {

Double3 VectorOnUnitCone(double phi, double theta) {
  // Height of vector in unit sphere
  double z = std::cos(theta);
  // Direct / shortest distance of point on sphere and z-axis
  double r = std::sqrt(1 - std::pow(z, 2.0));
  Double3 cone_vector_unit_sphere{r * std::cos(phi), r * std::sin(phi), z};
  return cone_vector_unit_sphere;
}

Double3 VectorOnConeAroundAxis(const Double3& axis, double phi, double theta) {
  // Create a vector in the unit sphere on cone (phi, theta) around (0,0,1)
  Double3 z_axis{0, 0, 1};
  Double3 cone_vector_unit_sphere = VectorOnUnitCone(phi, theta);

  // Rotate the vector form the unit sphere into the correct
  Double3 rotation_axis = Math::CrossProduct(z_axis, axis);
  double rotation_angle = std::acos(z_axis * axis);
  Double3 result = Math::RotAroundAxis(cone_vector_unit_sphere, rotation_angle,
                                       rotation_axis);
  return result;
}

}  // namespace bdm
