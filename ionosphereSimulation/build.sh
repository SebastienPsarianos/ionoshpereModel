#!/bin/bash
set -e

BUILD_SPACK=false

for arg in "$@"; do
    case "$arg" in
        --spack)
            BUILD_SPACK=true
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: $0 [--spack]"
            exit 1
            ;;
    esac
done

if [ -f "$HOME/spack/share/spack/setup-env.sh" ]; then
    . "$HOME/spack/share/spack/setup-env.sh"
else
    echo "Error: Could not find Spack setup script at ~/spack/share/spack/setup-env.sh"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

if [[ "$(uname)" == "Darwin" ]]; then
    SPACK_ENV_DIR="$SCRIPT_DIR/OSX"
else
    SPACK_ENV_DIR="$SCRIPT_DIR/linux"
fi

spack env activate "$SPACK_ENV_DIR"

if [ "$BUILD_SPACK" = true ]; then
    echo "Rebuilding Spack packages..."
    spack concretize -f --deprecated
    spack install --deprecated
fi

VIEW_PATH="$SPACK_ENV/.spack-env/view"
mkdir -p "$SCRIPT_DIR/build"
cd "$SCRIPT_DIR/build"
cmake .. \
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_PREFIX_PATH="$VIEW_PATH" \
    -DCMAKE_BUILD_TYPE=Release

make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)

echo "Build complete!"
