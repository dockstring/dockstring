#!/usr/bin/env bash

docker run -it --rm \
    --mount "type=bind,source=$(pwd),target=/workspaces/dockstring" \
    docker.io/library/dockstring:latest \
    bash