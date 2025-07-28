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


BINARY="${GPU_DIR}/bin/autodock_gpu_${NUMWI,,}wi"

### Step 1: Compile AutoDock-GPU
if [ ! -f "$BINARY" ]; then
    echo "[INFO] AutoDock-GPU binary not found. Compiling for $DEVICE..."

    if ! command -v nvcc &> /dev/null; then
        echo "[ERROR] CUDA toolkit not found. Please install it and try again."
        exit 1
    fi

    export GPU_INCLUDE_PATH="/usr/local/cuda/include"
    export GPU_LIBRARY_PATH="/usr/local/cuda/lib64"

    echo "[INFO] Running make DEVICE=$DEVICE NUMWI=$NUMWI..."
    make DEVICE=$DEVICE NUMWI=$NUMWI

    if [ -f "$BINARY" ]; then
        echo "[INFO] AutoDock-GPU compilation successful."
    else
        echo "[ERROR] AutoDock-GPU compilation failed."
        exit 2
    fi
else
    echo "[INFO] AutoDock-GPU binary already compiled."
fi

### Step 2: Compile AutoGrid
echo "[INFO] Checking AutoGrid..."

AUTOGRID_DIR="${GPU_DIR}/autogrid"
AUTOGRID_BINARY="${AUTOGRID_DIR}/autogrid4"

autogrid_needs_compile=false

# Check if binary exists and is runnable
if [ -x "$AUTOGRID_BINARY" ]; then
    if ! "$AUTOGRID_BINARY" --help &> /dev/null; then
        echo "[WARN] AutoGrid binary is not responding. Will recompile."
        autogrid_needs_compile=true
    else
        echo "[INFO] AutoGrid binary is functional."
    fi
else
    echo "[INFO] AutoGrid binary not found or not executable."
    autogrid_needs_compile=true
fi

# Compile if necessary
if [ "$autogrid_needs_compile" = true ]; then
    if [ ! -d "$AUTOGRID_DIR" ]; then
        echo "[ERROR] AutoGrid source directory not found at $AUTOGRID_DIR"
        exit 3
    fi

    cd "$AUTOGRID_DIR" || exit 4
    echo "[INFO] Compiling AutoGrid..."
    echo "[INFO] Preparing AutoGrid build environment..."
    if [ ! -f configure ]; then
        echo "[INFO] Running autoreconf to generate configure script..."
        autoreconf -i
    fi

    echo "[INFO] Running ./configure..."
    if [ ! -f configure ]; then
        echo "[INFO] No configure script found, running autoreconf..."
        autoreconf -i
    fi
    ./configure
    echo "[INFO] Cleaning previous build (if any)..."
    make clean || true
    echo "[INFO] Running make..."
    make -j$(nproc)

    if [ -x "$AUTOGRID_BINARY" ]; then
        echo "[INFO] AutoGrid compilation successful."
    else
        echo "[ERROR] AutoGrid compilation failed or binary not executable."
        exit 5
    fi

    cd "$GPU_DIR" || exit 1
fi

echo "All components are compiled and ready."
exit 0