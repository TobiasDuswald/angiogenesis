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

# This script is used to create a movie from a set of *.png files stored in a
# directory. The movie is saved in the same directory.
# usage: ./create_movie.sh <directory> <movie_name> <framerate>

set -e

# Get the director of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Default values for command line arguments
MOVIE_NAME="movie.mp4"
FRAMERATE="30"

# Load the command line arguments
DIRECTORY=$1
if [ -n "$2" ]; then
    MOVIE_NAME=$2
fi
if [ -n "$3" ]; then
    FRAMERATE=$3
fi

# Get cwd
CWD=$(pwd)

# Load the colors
source scripts/util.sh

# Check if the directory exists
if [ ! -d "$DIRECTORY" ]; then
    echo -e "${RED}Directory $DIRECTORY does not exist. ${NC}"
    exit 1
fi
cd $DIRECTORY

# Check if ffmpeg is installed
if ! command -v ffmpeg &> /dev/null
then
    echo -e "${RED}ffmpeg could not be found. Please install ffmpeg.${NC}"
    exit 1
fi

# Count the number of png files, if there are no png files, exit
NUM_PNG=$(ls -1 *.png 2>/dev/null | wc -l)
if [ $NUM_PNG == 0 ]; then
    echo -e "${RED}No png files found in $DIRECTORY. ${NC}"
    exit 1
else
    echo -e "${GREEN}<bash>${NC} Found $NUM_PNG png files in $DIRECTORY."
fi

# Create the movie if MOVIE_NAME does not exist
if [ -f "$MOVIE_NAME" ]; then
    echo -e "${RED}Movie $MOVIE_NAME already exists. Skipping it. ${NC}"
else
    echo -e "${GREEN}<bash>${NC} Creating movie $MOVIE_NAME with $FRAMERATE fps."
    ffmpeg -framerate $FRAMERATE -pattern_type glob -i '*.png' \
    -c:v libx264 -pix_fmt yuv420p $MOVIE_NAME
    echo -e "${GREEN}<bash>${NC} Done. Created movie $MOVIE_NAME."
fi

# Go back to the original directory
cd $CWD
