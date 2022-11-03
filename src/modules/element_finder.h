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

#ifndef ELEMENT_FINDER_H_
#define ELEMENT_FINDER_H_

#include <limits>
#include "modules/element_container.h"
#include "unibn_octree.h"

namespace bdm {

/// The TipCellFinder is a Octree based search engine to identify the closest
/// center of an element relative to a given vector.
class TipCellFinder {
 public:
  /// Constructing a TipCellFinder requires a mfem::Mesh. During this call, the
  /// Octree is build for fast searches.
  TipCellFinder();

  ~TipCellFinder() { delete octree_; };

  /// Calls back to FindClosestTipCell(mfem::Vector&).
  int FindClosestTipCell(const Real3& x) const;

  /// Get the center coordinates of a given element labeled by element_id.
  Real3 GetTipCellCenter(int element_id) const;

  /// Returns if there is at least one TipCell in the ball of radius r around
  /// the point x.
  bool IsTipCellInBall(const Real3& x, double r) const;

  size_t GetNumberOfTipCells() const { return tip_cell_container_.size(); }

  // Update the Octree
  void Update();

 private:
  /// Octree for spatial searches
  unibn::Octree<Real3, TipCellContainer>* octree_ = nullptr;

  /// Wraps the access to a mfem::Mesh for octree search
  TipCellContainer tip_cell_container_;

  /// Flag to avoid unnecessary updates of the container
  bool update_container_ = false;
};

}  // namespace bdm

#endif  // ELEMENT_FINDER_H_
