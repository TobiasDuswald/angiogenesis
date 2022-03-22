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

#include "mechanical_forces.h"
#include <math.h>
#include <stdexcept>

namespace bdm {

double MechanicalInteractionForce::EuclideanDistance(const Double3& x,
                                                     const Double3& y) const {
  return (x - y).Norm();
}

void MechanicalInteractionForce::CalculateAdhesiveForce(
    const Double3& direction, double distance, double sum_of_action_radi,
    Double3& result) const {
  if (distance > 0.0 && distance <= sum_of_action_radi) {
    double factor = -adhesion_scale_parameter_ *
                    std::pow(distance / sum_of_action_radi - 1, 2);
    result = direction;
    result *= factor;
  } else {
    result = {0.0, 0.0, 0.0};
  }
}

void MechanicalInteractionForce::CalculateRepulsiveForce(
    const Double3& direction, double distance, double sum_of_nuclear_radi,
    double sum_of_radi, Double3& result) const {
  double factor{0.0};
  if (distance > 0 && distance < sum_of_nuclear_radi) {
    factor = sum_of_nuclear_radi * distance / std::pow(sum_of_radi, 2) -
             2 * distance / sum_of_radi + 1;
  } else if (distance > sum_of_nuclear_radi && distance < sum_of_radi) {
    factor = std::pow(distance, 2) / std::pow(sum_of_radi, 2) -
             2 * distance / sum_of_radi + 1;
  } else {
    ;  // do nothing, keep factor as 0
  }
  factor *= repulsive_scale_parameter_;
  result = direction;
  result *= factor;
}

inline Double4 MechanicalInteractionForce::ConvertToDouble4(
    const Double3& x) const {
  Double4 y{x[0], x[1], x[2], 0};
  return y;
}

Double4 MechanicalInteractionForce::Calculate(const Agent* lhs,
                                              const Agent* rhs) const {
  // We cannot compute a interaction force between a Cell and itself
  if (lhs == rhs) {
    throw std::invalid_argument(
        "<MechanicalInteractionForce::Calculate>:"
        "cannot compute forces between an object and itself.");
  }
  // Cast Agents to TumorCells
  const TumorCell* lhs_tumor_cell = dynamic_cast<const TumorCell*>(lhs);
  const TumorCell* rhs_tumor_cell = dynamic_cast<const TumorCell*>(rhs);
  if (!lhs_tumor_cell or !rhs_tumor_cell) {
    // If neighbor is not a TumorCell, do not't calculate a force.
    return {0, 0, 0, 0};
  }

  // Initialize arrays for results
  Double3 result{0.0, 0.0, 0.0};
  Double3 adhesive_force{0.0, 0.0, 0.0};
  Double3 repulsive_force{0.0, 0.0, 0.0};

  // Get radius of lhs and rhs TumorCells
  const double radius_lhs = lhs_tumor_cell->GetRadius();
  const double radius_rhs = rhs_tumor_cell->GetRadius();
  const double nuclear_radius_lhs = lhs_tumor_cell->GetNuclearRadius();
  const double nuclear_radius_rhs = rhs_tumor_cell->GetNuclearRadius();
  const double action_radius_lhs = lhs_tumor_cell->GetActionRadius();
  const double action_radius_rhs = rhs_tumor_cell->GetActionRadius();

  // Compute distances and the vector direction of the connection. lhs - rhs is
  // necessary for compatibility with the sign of the forces.
  Double3 distance_vector =
      lhs_tumor_cell->GetPosition() - rhs_tumor_cell->GetPosition();
  double center_distance = distance_vector.Norm();
  distance_vector.Normalize();

  // Compute the forces
  CalculateAdhesiveForce(distance_vector, center_distance,
                         action_radius_lhs + action_radius_rhs, adhesive_force);
  CalculateRepulsiveForce(distance_vector, center_distance,
                          nuclear_radius_lhs + nuclear_radius_rhs,
                          radius_lhs + radius_rhs, repulsive_force);

  // Add up forces
  result = adhesive_force + repulsive_force;
  return ConvertToDouble4(result);
}

MechanicalInteractionForce* MechanicalInteractionForce::NewCopy() const {
  return new MechanicalInteractionForce(adhesion_scale_parameter_,
                                        repulsive_scale_parameter_);
}

}  // namespace bdm
