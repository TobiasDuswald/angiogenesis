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

#ifndef TIP_CELL_CONTAINER_H_
#define TIP_CELL_CONTAINER_H_
#include <vector>
#include "biodynamo.h"

namespace bdm {

class TipCellContainer {
 private:
  std::vector<Real3> element_center_points_;

 public:
  TipCellContainer();

  /// Copies the element center points of the mesh to element_center_points via
  /// an intermediate vector. Copies in serial.
  void Update();

  /// Returns the number of tip cells.
  size_t size() const;

  /// Retuns the center coordinate of a element in mesh_ labeled by idx.
  const Real3& operator[](size_t idx) const;
};

}  // namespace bdm

#endif  // TIP_CELL_CONTAINER_H_
