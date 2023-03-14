#!/bin/bash

# This script is used to run the full scale simulation.
# Usage: ./scripts/experiments/run_full_scale_model.sh

# Exit on error
set -e

# Get the director of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Get paths, compile BioDynaMo and the simulation, reset changes to bdm.json
source $DIR/compile.sh kFullScaleModel

# Change to the simulation directory
echo -e "${GREEN}<bash>${NC} Change to directory: $BDM_SCRIPT_DIR"
cd $BDM_SCRIPT_DIR

# Define parameters for the simulation
echo -e "${GREEN}<bash>${NC} Define parameters for the simulation"
if [ -f bdm.json ]; then
  mv bdm.json bdm.json.bak
fi
cp $BDM_SCRIPT_DIR/scripts/experiments/parameters/full_scale.json bdm.json

# Run the simulation (ToDo: Modify to meet experiment requirements)
echo -e "${GREEN}<bash>${NC} Run the simulation"
bdm run

# Reset changes to bdm.json
echo -e "${GREEN}<bash>${NC} Reset changes to bdm.json"
mv bdm.json.bak bdm.json

# Reset changes to src/angiogenesis_simulation.cc
echo -e "${GREEN}<bash>${NC} Reset changes to src/angiogenesis_simulation.cc"
git checkout $BDM_SCRIPT_DIR/src/angiogenesis_simulation.cc
