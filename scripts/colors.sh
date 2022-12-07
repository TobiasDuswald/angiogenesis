#!/bin/bash

# This script is used to define a list of colors that can be used in bash 
# scripts. 
# Usage: source scripts/colors.sh

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
