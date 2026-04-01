#!/bin/bash

# Navigate to the miracl_core/cpp directory
cd "$(dirname "$0")/miracl_core/cpp" || exit 1

# Automatically input the curve number and run the configuration script
echo -e "31\n0\n" | python3 config64.py