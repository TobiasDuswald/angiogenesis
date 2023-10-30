#!/bin/bash

# -----------------------------------------------------------------------------
#
# Copyright (C) 2022 CERN, TUM, and UT Austin. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# See the LICENSE file distributed with this work for details.
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# -----------------------------------------------------------------------------

# This script is used to run the vessel to center growth simulation.
# Usage: ./scripts/experiments/run_vessel_to_center.sh

# Exit on error
set -e

# Get the director of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Get paths, compile BioDynaMo and the simulation, reset changes to bdm.json
source $DIR/compile.sh kVesselsToCenter

# Change to the simulation directory
echo -e "${GREEN}<bash>${NC} Change to directory: $BDM_SCRIPT_DIR"
cd $BDM_SCRIPT_DIR

# Run the simulation (ToDo: Modify to meet experiment requirements)
echo -e "${GREEN}<bash>${NC} Run the simulation"
bdm run
