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

# This script is used to define a list of colors that can be used in bash 
# scripts. 
# Usage: source scripts/util.sh

# Define a list of colors
BLACK='\033[0;30m'
DARKGRAY='\033[1;30m'
RED='\033[0;31m'
LIGHTRED='\033[1;31m'
GREEN='\033[0;32m'
LIGHTGREEN='\033[1;32m'
ORANGE='\033[0;33m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
LIGHTBLUE='\033[1;34m'
PURPLE='\033[0;35m'
LIGHTPURPLE='\033[1;35m'
CYAN='\033[0;36m'
LIGHTCYAN='\033[1;36m'
LIGHTGRAY='\033[0;37m'
WHITE='\033[1;37m'
NC='\033[0m' # No Color

echo -e "${PURPLE}Co${BLUE}lo${LIGHTCYAN}rs ${GREEN}lo${YELLOW}ad${RED}ed${NC}"

# This function wraps the sed command to replace a string in a file. The
# function evaluates if we are running on a Mac or Linux and uses the
# appropriate sed command. The function takes three arguments:
# 1. The file to be modified
# 2. The string to be replaced
# 3. The string to replace the first string
function replace() {
    # Verify that the function is called with the correct number of arguments
    if [[ $# -ne 3 ]]; then
        echo -e "${RED}Error: replace() takes 3 arguments${NC}"
        echo -e "${RED}Usage: replace <file> <string to replace> <replacement>${NC}"
        return 1
    fi
    # Verify that file exists
    if [[ ! -f "$1" ]]; then
        echo -e "${RED}Error: File <$1> does not exist${NC}"
        return 1
    fi
    # Check if <string to replace> is in <file>
    if ! grep -q "$2" "$1"; then
        echo -e "${RED}Error: <$2> not found in $1${NC}"
        return 1htop
    fi
    # Call sed with the appropriate arguments
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "s/$2/$3/g" $1
    else
        sed -i "s/$2/$3/g" $1
    fi
    # Check if the replacement was successful
    if ! grep -q "$3" "$1"; then
        echo -e "${RED}Error: Replacement failed${NC}"
        return 1
    else
        echo -e "${GREEN}Success: <$2> replaced with <$3>${NC}"
    fi
    return 0
}
