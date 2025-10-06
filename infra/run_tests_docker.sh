#!/usr/bin/env bash

docker run -it --rm \
    --mount "type=bind,source=$(pwd),target=/code" \
    dockstring:latest \
    bash -c 'mamba activate dockstring && pip install pytest -e /code && pytest /code'
    
