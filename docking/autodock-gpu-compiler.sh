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

    # --------- ADDED: auto-detect TARGETS (keeps manual override if you set TARGETS/TARGET_ARCH) ---------
    detect_nvcc_support_set() {
      local major minor
      read major minor < <(nvcc --version | sed -n 's/.*release \([0-9]\+\)\.\([0-9]\+\).*/\1 \2/p')
      if [[ "$major" -ge 13 ]]; then
        echo "89 100 101 120"
      elif [[ "$major" -eq 12 && "$minor" -ge 8 ]]; then
        echo "89 100 101 120"
      elif [[ "$major" -eq 12 && "$minor" -ge 6 ]]; then
        echo "89 100 101"
      elif [[ "$major" -eq 12 ]]; then
        echo "80 86 89"
      else
        echo "50 52 53 60 61 62 70 72 75 80 86 89"
      fi
    }
    detect_gpu_targets() {
      if ! command -v nvidia-smi >/dev/null 2>&1; then
        echo ""
        return
      fi
      nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null \
        | awk '{gsub(/\./,""); print}' | sort -u | tr '\n' ' ' | sed 's/[[:space:]]*$//'
    }
    # Only set defaults if user hasn’t provided TARGETS/TARGET_ARCH already
    if [[ -z "${TARGETS:-}" || -z "${TARGET_ARCH:-}" ]]; then
      _detected="$(detect_gpu_targets)"       # e.g. "89 100 120"
      _supported="$(detect_nvcc_support_set)" # from nvcc version
      _choose=""
      if [[ -n "$_detected" ]]; then
        for t in $_detected; do
          if grep -qw "$t" <<<"$_supported"; then _choose+=" $t"; fi
        done
      fi
      if [[ -z "$_choose" ]]; then _choose="$_supported"; fi
      TARGETS_DEFAULT="$(echo "$_choose" | xargs)"  # trim spaces
      # prefer the highest for TARGET_ARCH (first item after sorting numerically desc)
      TARGET_ARCH_DEFAULT="sm_$(for t in $TARGETS_DEFAULT; do echo $t; done | sort -nr | head -n1)"
      # export only if not preset
      TARGETS="${TARGETS:-$TARGETS_DEFAULT}"
      TARGET_ARCH="${TARGET_ARCH:-$TARGET_ARCH_DEFAULT}"
    fi
    echo "[INFO] Auto-detected TARGETS: ${TARGETS} (TARGET_ARCH=${TARGET_ARCH})"
    # --------- end ADDED ---------

    echo "[INFO] Running make DEVICE=$DEVICE NUMWI=$NUMWI..."
    mkdir -p "${GPU_DIR}/bin"
    make DEVICE=$DEVICE NUMWI=$NUMWI TARGET_ARCH="${TARGET_ARCH}" TARGETS="${TARGETS}" NVCC="nvcc" CC=gcc CXX=g++   ### USE CUDATOOLKIT 12.8, CUDA TOOLKIT 13 IS NOT SUPPORTED YET
    #####PLEASE ADJUST THE TARGET_ARCH AND TARGETS AS NEEDED FOR YOUR GPU
    ####Guidelines for choosing the right target architecture:
    # - sm_80 for Volta GPUs (e.g., Tesla V100)
    # - sm_86 for Ampere GPUs (e.g., A100, RTX 3000 series)
    # - sm_89 for Ada Lovelace GPUs (e.g., RTX 4000 series)
    # - sm_90 for Hopper GPUs (e.g., H100 variants)
    # - sm_100/101 for Blackwell GB-series (datacenter)
    # - sm_120 for consumer Blackwell (e.g., some RTX 50xx)
    ### you got the idea, just adjust the TARGET_ARCH and TARGETS variables.
    TODO=$(grep -E 'TARGET_ARCH|TARGETS' Makefile.Cuda)
    echo "[INFO] Makefile.Cuda settings: $TODO"

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
