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

# This script is used to extract the results of the application. It duplicates
# the <input_dir> and stores the results in <output_dir>. It copies all files
# except for the file types that are not needed (*.pvsm, *.vti, *.vtu, *.pvti,
# *.pvtu). It then compresses the files in the output directory, and deletes the
# uncompressed files (<output_dir>).
# usage: ./extract_results.sh <input_dir> <output_dir>

set -e

# Get the path of this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source $DIR/util.sh

# Check if the number of arguments is correct
if [ $# -ne 2 ]; then
    echo -e "${RED}<bash>${NC} Usage: ./extract_results.sh <input_dir> <output_dir>"
    exit 1
fi

# Define the input directory
INPUT_DIR=$1

# Define the output directory
OUTPUT_DIR=$2

# Check if the input directory exists
if [ ! -d $INPUT_DIR ]; then
    echo -e "${RED}<bash>${NC} Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Check if the output directory does exist, if it does, exit the script
if [ -d $OUTPUT_DIR ]; then
    echo -e "${RED}<bash>${NC} Output directory already exists: $OUTPUT_DIR"
    exit 1
fi

# Use rsync to copy the files from the input directory to the output directory
# but exclude the file types that are not needed (*.pvsm, *.vti, *.vtu, *.pvti,
# *.pvtu)
rsync -av --exclude="*.pvsm" --exclude="*.vti" --exclude="*.vtu" --exclude="*.pvti" --exclude="*.pvtu" $INPUT_DIR/ $OUTPUT_DIR

# Compres the output directory
tar -zcvf $OUTPUT_DIR.tar.gz $OUTPUT_DIR

# Remove the output directory
rm -rf $OUTPUT_DIR

# Print the path to the compressed output directory
echo -e "${GREEN}<bash>${NC} Compressed output directory: $OUTPUT_DIR.tar.gz"
