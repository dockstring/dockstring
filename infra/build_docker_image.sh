#!/usr/bin/env bash

# exit if any command fails
set -e 
set -o pipefail

BASE_NAME="dockstring"
TAG=${BASE_NAME}":"$(date +"%Y%m%dT%H%M%S")
FILE="infra/dockstring.dockerfile"

echo "Building Docker image with tag \"${TAG}\"" 

docker build \
    --progress=plain \
    --tag=${TAG} \
    --file=${FILE} \
    . \
    2>&1 | tee -a "build_${TAG}.log"

docker tag ${TAG} ${BASE_NAME}":latest"
