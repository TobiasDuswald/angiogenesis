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

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "biodynamo.h"

namespace bdm {

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

/// Returns the unit vector defined by phi and theta in spherical coordinates.
/// Helper function for VectorOnConeAroundAxis.
Double3 VectorOnUnitCone(double phi, double theta);

/// Returns a unit vector that lies on the cone defined by (phi, theta) around
/// axis. Intended for behaviors where phi is random.
Double3 VectorOnConeAroundAxis(const Double3& axis, double phi, double theta);

/// Takes three vectors and fills the latter two such that they are orthogonal
/// and normalized.
void GetOrthogonalSystem(const Double3& a, Double3& b, Double3& c);

// ---------------------------------------------------------------------------
// Sorting
// Credit:
// https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-\
// the-same-way-with-criteria-that-uses-only-one-of
// Usage example in test/test_suit_vector_operations.cc
// ---------------------------------------------------------------------------

/// Get a permutation that sorts the vector.
template <typename T, typename Compare>
std::vector<std::size_t> GetSortPermutation(const std::vector<T>& vec,
                                            const Compare& compare) {
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) {
    return compare(vec[i], vec[j]);
  });
  return p;
};

// Apply permutation to vector (in place).
template <typename T>
void ApplyPermutationInPlace(std::vector<T>& vec,
                             const std::vector<std::size_t>& p) {
  std::vector<bool> done(vec.size());
  for (std::size_t i = 0; i < vec.size(); ++i) {
    if (done[i]) {
      continue;
    }
    done[i] = true;
    std::size_t prev_j = i;
    std::size_t j = p[i];
    while (i != j) {
      std::swap(vec[prev_j], vec[j]);
      done[j] = true;
      prev_j = j;
      j = p[j];
    }
  }
};

}  // namespace bdm

#endif  // VECTOR_OPERATIONS_H_
