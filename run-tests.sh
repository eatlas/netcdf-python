#!/usr/bin/env bash

# Run the Python tests using unittest.

# Identify the directory containing this executable file, as that is where the Dockerfile will be found.
SCRIPT=$(readlink -f "$0")
SCRIPT_PATH=$(dirname "${SCRIPT}")

python -m unittest discover -s ${SCRIPT_PATH}/test -t ${SCRIPT_PATH}
