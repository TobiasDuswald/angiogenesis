#!/bin/bash

# Exit on failure
set -e

# Define colors: Orange, green, normal
ORANGE='\033[0;33m'
GREEN='\033[0;32m'
NC='\033[0m'

# Get the name of the script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Zip the folder results and add a timestamp to the filename
echo -e "${ORANGE}Zipping results folder .. ${NC}"
ZIPNAME=results-$(date +%Y-%m-%d-%H-%M-%S).zip
zip -r $ZIPNAME results
echo -e "${GREEN}Zipping results folder .. done${NC}"

# Define Cernbox folder
CERNBOX_FOLDER=~/Cernbox/BDM/Angiogenesis

# Create the CERNBOX folder recursively if it does not exist
mkdir -p $CERNBOX_FOLDER

# Copy the zip file to the backup folder ~/Cernbox
echo -e "${ORANGE}Copy $ZIPNAME to $CERNBOX_FOLDER ${NC}"
cp $DIR/$ZIPNAME $CERNBOX_FOLDER
echo -e "${GREEN}Copy $ZIPNAME to $CERNBOX_FOLDER .. done${NC}"

# Remove the zip file
echo -e "${ORANGE}Cleanup .. ${NC}"
rm $DIR/$ZIPNAME
echo -e "${GREEN}Cleanup .. done ${NC}"

# Exit the script
echo -e "\n${GREEN}Backup finished successfully."
exit 0
