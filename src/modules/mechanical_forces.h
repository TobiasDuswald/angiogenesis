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

/// ----------------------------------------------------------------------------
/// In this file, we define a custom interaction force between TumorCells. The
/// force is defined as in Rocha et al 2018 / Lima et al 2021 and takes into
/// account repulsion and adhesion.
/// ----------------------------------------------------------------------------

#ifndef MECHANICAL_FORCES_H_
#define MECHANICAL_FORCES_H_

#include "biodynamo.h"
#include "modules/tumor_cell.h"

namespace bdm {

// Implementation of the specific force module representing the force suggested
// in the work from Lima et al.
class MechanicalInteractionForce : public InteractionForce {
 private:
  // Compute the Euclidean (L2) distance between two points x and y.
  double EuclideanDistance(const Double3& x, const Double3& y) const;

  // Computes the adhesive force vector `result` of two cells at positions x and
  // y.
  // @param[in] direction: direction of the force: (x-y)/||x-y||
  // @param[in] distance: distance between the two cells: ||x-y||
  // @param[in] sum_of_action_radi: the sum of the action radi of cell 1 and 2
  // @param[in/out] result: vector of adhesive force
  void CalculateAdhesiveForce(const Double3& direction, double distance,
                              double sum_of_action_radi, Double3& result) const;

  // Computes the repulsive force vector `result` of two cells at positions x
  // and y.
  // @param[in] direction: direction of the force: (x-y)/||x-y||
  // @param[in] distance: distance between the two cells: ||x-y||
  // @param[in] sum_of_action_radi: the sum of the action radi of cell 1 and 2
  // @param[in/out] result: vector of repulsive force
  void CalculateRepulsiveForce(const Double3& direction, double distance,
                               double sum_of_nuclear_radi, double sum_of_radi,
                               Double3& result) const;

  // Converts a Double3 x to Double4 by returning {x0,x1,x2,0.0}
  Double4 ConvertToDouble4(const Double3& x) const;

  // Numeric parameter c_{cca} for force, unit [\mu m / min]
  const double adhesion_scale_parameter_;
  // Numeric parameter c_{ccr}for force, unit [\mu m / min]
  const double repulsive_scale_parameter_;

 public:
  MechanicalInteractionForce() = delete;
  MechanicalInteractionForce(double adhesion_scale_parameter,
                             double repulsive_scale_parameter)
      : adhesion_scale_parameter_(adhesion_scale_parameter),
        repulsive_scale_parameter_(repulsive_scale_parameter){};
  virtual ~MechanicalInteractionForce() {}

  // Calculate the interaction force between two agents.
  Double4 Calculate(const Agent* lhs, const Agent* rhs) const override;

  MechanicalInteractionForce* NewCopy() const override;
};

}  // namespace bdm

#endif  // MECHANICAL_FORCES_H_
