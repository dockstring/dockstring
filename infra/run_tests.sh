#!/usr/bin/env bash

docker run -it --rm \
    --mount "type=bind,source=$(pwd),target=/dockstring" \
    docker.io/library/dockstring:latest \
    bash -c "conda activate dockstring && pip install /dockstring pytest && pytest /dockstring --collect-only"
    