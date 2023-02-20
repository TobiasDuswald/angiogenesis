#!/bin/bash

# This script is used to run the full scale simulation.
# Usage: ./scripts/experiments/run_full_scale_model.sh

# Exit on error
set -e

# Get the director of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Get paths, compile BioDynaMo and the simulation, reset changes to bdm.json
source $DIR/compile.sh

# Change to the simulation directory
echo -e "${GREEN}<bash>${NC} Change to directory: $BDM_SCRIPT_DIR"
cd $BDM_SCRIPT_DIR

# Run the simulation (ToDo: Modify to meet experiment requirements)
echo -e "${GREEN}<bash>${NC} Run the simulation"
bdm run
