# syntax=docker/dockerfile:1
ARG UBUNTU_VERSION=bionic-20230530  # for repro-builds
FROM ubuntu:${UBUNTU_VERSION}

RUN apt-get update && apt-get install -y \
    wget \
    # git is being used in the package
    git \
    && rm -rf /var/lib/apt/lists/*

SHELL [ "/bin/bash", "-c" ]

# mamba installation
ENV MAMBA_DIR=/opt/miniforge3
ENV PATH=${PATH}:${MAMBA_DIR}/bin
ENV MAMBA_ROOT_PREFIX=${MAMBA_DIR}

RUN wget --no-hsts --quiet --output-document=miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash miniforge.sh -b -p ${MAMBA_DIR} \
    && rm miniforge.sh \
    && unlink ${MAMBA_DIR}/bin/python3.1 \
    && ${MAMBA_DIR}/bin/mamba update "mamba>2.0" \
    && echo "eval $(mamba shell hook --shell bash)" >> /etc/profile.d/source_mamba.sh \
    # for interactive shells:
    && echo "source /etc/profile.d/source_mamba.sh" >> /etc/bash.bashrc

# for non-interactive, not login shells:
# https://www.solipsys.co.uk/images/BashStartupFiles1.png
ENV BASH_ENV="/etc/profile.d/source_mamba.sh"

# create environment
COPY ./infra/repro-environment.yml ./environment.yml
RUN mamba env create --file environment.yml \
    && mamba clean --all --yes \
    && rm environment.yml
