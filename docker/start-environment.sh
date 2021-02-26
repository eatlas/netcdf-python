#!/usr/bin/env bash

# Start a Docker container for development.
#
# Usage:
#
#    start-environment.sh <data directory on host machine>

SCRIPT=$(readlink -f "$0")
SCRIPT_PATH=$(dirname "${SCRIPT}")
ROOT_PATH=$(readlink -f "${SCRIPT_PATH}/..")

# Verify the data directory has been specified.
if [ -z "$1" ]
  then
    echo "Usage: start-environment.sh <data directory on host machine>"
    exit 1
fi
if [ ! -d "$1" ]
  then
    echo "Specified data directory \"$1\" does not exist."
    exit 1
fi

# Pull the latest Python 3 Docker Image.
docker pull python:3

# Build the Docker image.
docker build -f ${SCRIPT_PATH}/Dockerfile -t python-netcdf:dev ${SCRIPT_PATH}

# Run the Docker image with relevant mount points.
docker run \
  -it \
  --rm \
  -u $(id -u):$(id -g) \
  -w /workdir \
  -v ${ROOT_PATH}:/workdir \
  python-netcdf:dev bash
