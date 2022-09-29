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

#ifndef COUNT_NEIGHBORS_H_
#define COUNT_NEIGHBORS_H_

#include "biodynamo.h"

namespace bdm {

template <class AgentToCount>
class CountNeighborsFunctor : public Functor<void, Agent*, double> {
 private:
  size_t num_neighbors_;
  double squared_distance_;

 public:
  CountNeighborsFunctor(double distance)
      : num_neighbors_(0), squared_distance_(std::pow(distance, 2.0)) {}

  // This is called once for each neighbor that is found
  void operator()(Agent* neighbor, double squared_distance) {
    if (dynamic_cast<AgentToCount>(neighbor) &&
        squared_distance < squared_distance_)
#pragma omp atomic
      num_neighbors_ += 1;
  }

  double GetNumNeighbors() { return num_neighbors_; }

  void Reset() { num_neighbors_ = 0; }
};

}  // namespace bdm

#endif  // COUNT_NEIGHBORS_H_
