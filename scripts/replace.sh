#!/bin/bash

# This function wraps the sed command to replace a string in a file. The
# function evaluates if we are running on a Mac or Linux and uses the
# appropriate sed command. The function takes three arguments:
# 1. The file to be modified
# 2. The string to be replaced
# 3. The string to replace the first string
function replace() {
    # Verify that the function is called with the correct number of arguments
    if [[ $# -ne 3 ]]; then
        echo "Usage: replace <file> <string to replace> <replacement string>"
        return 1
    fi
    # Check if <string to replace> is in <file>
    if ! grep -q "$2" "$1"; then
        echo "String '$2' not found in file '$1'"
        return 1
    fi
    # Call sed with the appropriate arguments
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "s/$2/$3/g" $1
    else
        sed -i "s/$2/$3/g" $1
    fi
    # Check if the replacement was successful
    if ! grep -q "$3" "$1"; then
        echo "Replacement failed"
        return 1
    fi
    return 0
}
