#!/bin/bash
# This script is used to trim a set of images. We use imagemagick to trim the
# images, e.g remove the background around the image of interest. Use this
# script to, for instance, prepare images for the document focusing on the 
# tumor cells only, rather than having large empty spaces around the tumor.
# If imagemagick is not installed, the script will exit.
# Usage: ./trim.sh INPUT_DIR

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
INPUT_DIR=$1

# If imagemagick is not installed, exit the script
if ! command -v convert &> /dev/null
then
    echo -e "${RED}imagemagick could not be found. Please install imagemagick.${NC}"
    exit 0
fi

# Get cwd
CWD=$(pwd)

# Change into INPUT_DIR
cd $INPUT_DIR

# Print the number of png files in the directory
echo -e "${GREEN}<bash>${NC} Number of png files: $(ls -1 *.png 2>/dev/null | wc -l)"

# Create a folder for the cropped images, if it exists, delete it
if [ -d "trimmed" ]; then
    rm -rf trimmed
fi
mkdir trimmed

# Trim each image in the current folder
for file in *.png
do
    echo -e "${GREEN}<bash>${NC} Trim $file"
    convert $file -trim trimmed/$file
done

# Print the number of trimmed png files in the directory
echo -e "${GREEN}<bash>${NC} Number of trimmed png files: $(ls -1 trimmed/*.png 2>/dev/null | wc -l)"

cd $CWD
