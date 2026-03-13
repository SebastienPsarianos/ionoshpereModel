#!/bin/bash
set -e

cd "$(cd "$(dirname "$0")" && pwd)"

if [ -f "$HOME/spack/share/spack/setup-env.sh" ]; then
    . "$HOME/spack/share/spack/setup-env.sh"
else
    echo "Error: Could not find Spack setup script at ~/spack/share/spack/setup-env.sh"
    exit 1
fi

if [[ "$(uname)" == "Darwin" ]]; then
    SPACK_ENV_DIR="OSX"
else
    SPACK_ENV_DIR="linux"
fi

spack env activate "$SPACK_ENV_DIR"

mpirun -np 8 ./build/IonosphereSolver "$@"
