#!/bin/bash
# Install Julia via conda-forge if not available
if ! command -v julia &>/dev/null; then
    echo "Installing Julia via conda-forge..."
    if ! command -v micromamba &>/dev/null; then
        curl -fsSL -o /tmp/micromamba "https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64"
        chmod +x /tmp/micromamba
    fi
    export MAMBA_ROOT_PREFIX=/tmp/mamba
    /tmp/micromamba create -n julia -c conda-forge julia -y 2>&1
    echo 'export MAMBA_ROOT_PREFIX=/tmp/mamba' >> ~/.bashrc
    echo 'eval "$(/tmp/micromamba shell hook -s bash)"' >> ~/.bashrc
    echo 'micromamba activate julia' >> ~/.bashrc
fi
