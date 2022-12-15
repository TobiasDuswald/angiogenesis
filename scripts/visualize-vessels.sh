#!/bin/bash

# This script is used to visualize the output of the application. We expect that
# the output of the application is stored in the output directory. This script
# will find all paraview state files (*.pvsm) in the output directory and run
# pysrc/paraview_tumor_rotation.py with the state file.

set -e

# Define parameters
BACKGROUND=1 # (0: regular, 1: transparent)
AXES=0 # (0: off, 1: on)
COLORBAR=0 # (0: off, 1: on)
OVERWRITE=1 # (0: off, 1: on)

# Get the director of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Source colors
. $DIR/colors.sh

# Check if BioDynaMo is sourced by checking if the environment variable 
# BDMSYS is set
if [ -z "$BDMSYS" ]; then
    echo "BioDynaMo is not sourced. Please source BioDynaMo before running this script."
    exit 1
fi

# Define the path to the pvpython 
PVPYTHON="$BDMSYS/third_party/paraview/bin/pvpython"

# Define the output directory
OUTPUT_DIR="$DIR/../output"

# Get cwd
CWD=$(pwd)

# Change into the parent directory of OUTPUT_DIR
cd $OUTPUT_DIR/..

# Find all paraview state files (*.pvsm) in the output directory
FILES=$(find $OUTPUT_DIR -name "*.pvsm")
for file in $FILES
do
    # Get a timestamp
    timestamp=$(date +%s)
    # Extract the absolute path of the file
    file=$(realpath $file)
    # Run paraview with the state file
    # paraview $file
    echo -e "${GREEN}<bash>${NC} State: $file"
    echo -e "${GREEN}<bash>${NC} Render vessels"
    $PVPYTHON $DIR/../pysrc/paraview_vegf_vessels.py $file $BACKGROUND $AXES $COLORBAR $OVERWRITE
    # Print elapsed time
    echo -e "${GREEN}<bash>${NC} Elapsed time: $(($(date +%s)-timestamp)) seconds"
done

cd $CWD

# If ffmpeg is not installed, exit the script
if ! command -v ffmpeg &> /dev/null
then
    echo -e "${RED}ffmpeg could not be found. Please install ffmpeg.${NC}"
    exit 0
fi

# Iterate over all files and call the script create_movie.sh with the directory
# of the state file as argument
for file in $FILES
do
    # Get the directory of the state file
    dir=$(dirname $file)
    # Call the script create_movie.sh with the directory of the state file as
    # argument
    $DIR/create_movie.sh $dir/vessel_vegf_bg${BACKGROUND}_cb${COLORBAR}_ax${AXES} vessel_vegf.mp4 15
    $DIR/create_movie.sh $dir/vessel_bg${BACKGROUND}_cb${COLORBAR}_ax${AXES} vessel.mp4 15
done
