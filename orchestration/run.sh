#!/bin/bash
set -e

cd "$(dirname "$0")"

if [ -f "$HOME/spack/share/spack/setup-env.sh" ]; then
    . "$HOME/spack/share/spack/setup-env.sh"
else
    echo "Error: Could not find Spack setup script at ~/spack/share/spack/setup-env.sh"
    exit 1
fi

if [[ "$(uname)" == "Darwin" ]]; then
    SPACK_ENV_DIR="../ionosphereSimulation"
else
    SPACK_ENV_DIR="../ionosphereSimulation/spack-linux"
fi

spack env activate "$SPACK_ENV_DIR"

python3 orchestrator.py "$@"
