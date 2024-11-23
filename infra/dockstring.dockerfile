# syntax=docker/dockerfile:1
FROM ubuntu:bionic-20230530

RUN apt-get update && apt-get install -y \
    wget \
    # git is being used in the package
    git \
    && rm -rf /var/lib/apt/lists/*

SHELL [ "/bin/bash", "-c" ]

# mamba/conda installation
ENV CONDA_DIR=/opt/miniforge3
RUN wget --no-hsts --quiet --output-document=miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    && bash miniforge.sh -b -p ${CONDA_DIR} \
    && rm miniforge.sh \
    && unlink ${CONDA_DIR}/bin/python3.1 \
    && echo "source ${CONDA_DIR}/etc/profile.d/conda.sh" >> /etc/profile.d/source_conda.sh \
    # for interactive shells:
    && echo "source /etc/profile.d/source_conda.sh" >> /etc/bash.bashrc

# for non-interactive, not login shells:
# https://www.solipsys.co.uk/images/BashStartupFiles1.png
ENV BASH_ENV="/etc/profile.d/source_conda.sh"

# create environment
COPY ./infra/repro-environment.yml ./environment.yml
RUN conda env create --file environment.yml \
    && conda clean --all --yes \
    && rm environment.yml