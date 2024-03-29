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

#include "transition_probabilities.h"
#include <cmath>

namespace bdm {

double ComputeProbability_Q_To_SG2(const double nutrients, const double tra,
                                   const double delta_t,
                                   const SimParam* sparam) {
  return P_Q_SG2_N(nutrients, delta_t, sparam) *
         P_Q_SG2_TRA(tra, delta_t, sparam);
};

double ComputeProbability_SG2_To_SG2(const double dox, const double delta_t,
                                     const SimParam* sparam) {
  return P_SG2_SG2_DOX(dox, delta_t, sparam);
};

double ComputeProbability_SG2_To_D(const double dox, const double delta_t,
                                   const SimParam* sparam) {
  return P_SG2_D_DOX(dox, delta_t, sparam);
};

double ComputeProbability_Q_To_D(const double nutrients, const double tra,
                                 const double dox, const double delta_t,
                                 const SimParam* sparam) {
  return P_Q_D_N(nutrients, delta_t, sparam) *
         (1 + sparam->zeta_Q_D_TRA * tra + sparam->zeta_Q_D_DOX * dox +
          sparam->zeta_Q_D_TRA_DOX * tra * dox);
};

double ComputeProbability_H_To_D(const double tra, const double dox,
                                 const double delta_t, const SimParam* sparam) {
  return sparam->base_rate_H_D * delta_t *
         (1 + sparam->zeta_H_D_TRA * tra + sparam->zeta_H_D_DOX * dox +
          sparam->zeta_H_D_TRA_DOX * tra * dox);
  ;
};

double P_Q_SG2_N(const double nutrients, const double delta_t,
                 const SimParam* sparam) {
  return LinearProbabilityIncreaseForConcentration(
      nutrients, sparam->threshold_Q_SG2_N, sparam->alpha_Q_SG2_N, delta_t);
};

double P_Q_SG2_TRA(const double tra, const double delta_t,
                   const SimParam* sparam) {
  return std::exp(-sparam->alpha_Q_SG2_TRA * tra);
};

double P_Q_D_N(const double nutrients, const double delta_t,
               const SimParam* sparam) {
  return SmoothHeavisideForConcentration(nutrients, sparam->threshold_Q_D_N,
                                         sparam->alpha_Q_D_N, sparam->k_Q_D_N,
                                         delta_t, sparam->gamma_Q_D_N);
};

double P_Q_D_TRA(const double tra, const double delta_t,
                 const SimParam* sparam) {
  return 1 + sparam->zeta_Q_D_TRA * tra;
};

double P_Q_D_DOX(const double dox, const double delta_t,
                 const SimParam* sparam) {
  return 1 + sparam->zeta_Q_D_DOX * dox;
};

double P_SG2_SG2_DOX(const double dox, const double delta_t,
                     const SimParam* sparam) {
  return LinearProbabilityIncreaseForConcentration(
      dox, sparam->threshold_SG2_SG2_DOX, sparam->alpha_SG2_SG2_DOX, delta_t);
};

double P_SG2_D_DOX(const double dox, const double delta_t,
                   const SimParam* sparam) {
  return LinearProbabilityIncreaseForConcentration(
      dox, sparam->threshold_SG2_D_DOX, sparam->alpha_SG2_D_DOX, delta_t);
};

double P_H_D_DOX(const double dox, const double delta_t,
                 const SimParam* sparam) {
  return 1 + sparam->zeta_H_D_DOX * dox;
};

double P_H_D_TRA(const double tra, const double delta_t,
                 const SimParam* sparam) {
  return 1 + sparam->zeta_H_D_TRA * tra;
};

}  // namespace bdm
