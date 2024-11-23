#!/usr/bin/env bash

docker run -it --rm \
    --mount "type=bind,source=$(pwd),target=/dockstring" \
    dockstring:latest \
    bash -c "conda activate dockstring && pip install /dockstring pytest && pytest /dockstring"
    