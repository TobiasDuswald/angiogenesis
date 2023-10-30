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

# This script wraps the trim.sh script to apply it to all subfolders of a base
# directory. The script will trim all images in the subfolders and save the
# trimmed images in a folder called trimmed in the subfolder.
# Usage: ./trim_all.sh INPUT_DIR

set -e

# Get the director of the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Source colors
. $DIR/util.sh

# Check if we receive exactly one argument and save it in the variable INPUT_DIR
if [ $# -ne 1 ]; then
    echo -e "${RED}Usage: $0 INPUT_DIR${NC}"
    exit 1
fi

# Get the input directory
INPUT_DIR=$1

# Get the absolute path of the input directory
INPUT_DIR=$(realpath $INPUT_DIR)

# Get all subfolders of the input directory that conatin the substring 
# "cb0_ax0" and save them in the variable SUBFOLDERS
SUBFOLDERS=$(find $INPUT_DIR -type d -name "*cb0_ax0*")

# Loop over all subfolders
for SUBFOLDER in $SUBFOLDERS
do
    # Print the subfolder
    echo -e "${GREEN}<bash>${NC} Subfolder: $SUBFOLDER"
    # Trim the images in the subfolder
    $DIR/trim.sh $SUBFOLDER
done
