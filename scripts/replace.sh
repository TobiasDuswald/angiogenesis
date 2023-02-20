#!/bin/bash

# This function wraps the sed command to replace a string in a file. The
# function evaluates if we are running on a Mac or Linux and uses the
# appropriate sed command. The function takes three arguments:
# 1. The file to be modified
# 2. The string to be replaced
# 3. The string to replace the first string
function replace() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "s/$2/$3/g" $1
    else
        sed -i "s/$2/$3/g" $1
    fi
}
