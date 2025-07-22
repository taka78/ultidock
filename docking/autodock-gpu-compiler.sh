#!/bin/bash

# Takes the AUTODOCK_GPU_DIR as argument or uses current directory
GPU_DIR="${1:-$(pwd)}"
cd "$GPU_DIR" || exit 1

DEVICE="${DEVICE:-CUDA}"       # Default device type
NUMWI="${NUMWI:-128}"           # Default work items

###autodock_gpu_<N>wi refers to the number of OpenCL/CUDA work-items per ligand 
#(like thread-group size):
#64wi = 64 threads per ligand
#128wi = 128 threads per ligand (more parallelism)
#BINARY="./bin/autodock_gpu_${NUMWI,,}wi"

# Check if binary exists
if [ ! -f "$BINARY" ]; then
    echo "[INFO] AutoDock-GPU binary not found. Compiling for $DEVICE..."

    if [ "$DEVICE" == "CUDA" ]; then
        if ! command -v nvcc &>/dev/null; then
            echo "[ERROR] CUDA toolkit not found. Please install it and try again."
            exit 1
        fi
        export GPU_INCLUDE_PATH="/usr/local/cuda/include"
        export GPU_LIBRARY_PATH="/usr/local/cuda/lib64"
    fi

    make DEVICE="$DEVICE" NUMWI="$NUMWI" -C "$GPU_DIR"

    if [ ! -f "$BINARY" ]; then
        echo "[ERROR] Compilation failed."
        exit 2
    fi
else
    echo "[INFO] AutoDock-GPU binary already exists at $BINARY"
fi

# Run it if you passed extra args
if [ "$#" -gt 1 ]; then
    shift
    "$BINARY" "$@"
fi
