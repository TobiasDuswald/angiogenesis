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

#ifndef TIP_CELL_FINDER_OPERATION_H_
#define TIP_CELL_FINDER_OPERATION_H_

#include "angiogenesis_simulation.h"
#include "biodynamo.h"

namespace bdm {

// Operation to update the tip cell finder.
class UpdateTipCellFinder : public StandaloneOperationImpl {
  BDM_OP_HEADER(UpdateTipCellFinder);

  // Apply operation
  void operator()() override {
    auto* sim = Simulation::GetActive();
    // WARNING: This is a bit of a hack. We need to cast the simulation to the
    // derived class AngiogenesisSimulation to access the tip cell finder.
    // The simulation is however not polymorphic, so we cannot use dynamic_cast.
    // We can however use static_cast, because we know that the simulation is
    // an AngiogenesisSimulation. This is only true for this specific example.
    auto* asim = static_cast<AngiogenesisSimulation*>(sim);
    asim->UpdateTipCellFinder();
  };
};

}  // namespace bdm

#endif  // TIP_CELL_FINDER_OPERATION_H_
