#!/bin/bash

# This script is used to compile BioDynaMo and the siulation. If you don't want
# to compile BioDynaMo, then you can set the environment variable
# BPE_NO_CLEANBUILD to 1. Usage:
# ./scripts/experiments/compile.sh
# or when including the script in another script:
# source ./scripts/experiments/compile.sh <ExperimentID>

# Exit on error
set -e

# Get the path to the current directory and to the project root directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
BDM_SCRIPT_DIR=$DIR/../..
BDM_SCRIPT_DIR=$(realpath $BDM_SCRIPT_DIR)
EXAMPLES_DIR=$BDM_SCRIPT_DIR/../..
EXAMPLES_DIR=$(realpath $EXAMPLES_DIR)
PROJECT_ROOT_DIR=$BDM_SCRIPT_DIR/../../..
PROJECT_ROOT_DIR=$(realpath $PROJECT_ROOT_DIR)

# Source some utility scripts
source $BDM_SCRIPT_DIR/scripts/util.sh

# Print previously defined variables
echo -e "${GREEN}<bash>${NC} BDM_SCRIPT_DIR: $BDM_SCRIPT_DIR"
echo -e "${GREEN}<bash>${NC} EXAMPLES_DIR: $EXAMPLES_DIR"
echo -e "${GREEN}<bash>${NC} PROJECT_ROOT_DIR: $PROJECT_ROOT_DIR"

# If the project is a git repository, then reset all chanages to bdm.json, else
# exit with an error
if [ -d "$BDM_SCRIPT_DIR/.git" ]; then
    echo -e "${GREEN}<bash>${NC} Reset changes to bdm.json"
    git checkout $BDM_SCRIPT_DIR/bdm.json
else
    echo -e "${RED}Error: The project is not a git repository.${NC}"
    exit 1
fi

# Source the main utility script
echo -e "${GREEN}<bash>${NC} Source ${EXAMPLES_DIR}/util/main.sh"
source $EXAMPLES_DIR/util/main.sh

# Change to BDM_SCRIPT_DIR
echo -e "${GREEN}<bash>${NC} Change to $BDM_SCRIPT_DIR"

# Replace "experiment = Experiment::kAvascularTumorSpheroid" with 
# "experiment = Experiment::<ExperimentID>" in src/angiogenesis_simulation.cc
echo -e "${GREEN}<bash>${NC} Set experiment to $1"
replace $BDM_SCRIPT_DIR/src/angiogenesis_simulation.cc "experiment = Experiment::kAvascularTumorSpheroid" "experiment = Experiment::$1"

# Build biodynamo and the simulation
echo -e "${GREEN}<bash>${NC} Build biodynamo and the simulation"
source $EXAMPLES_DIR/util/default-compile-script.sh "-DCMAKE_BUILD_TYPE=Release -Djemalloc=off" "-DCMAKE_BUILD_TYPE=Release -Djemalloc=off"
