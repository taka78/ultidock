#!/bin/bash

set -euo pipefail

# Takes the AUTODOCK_GPU_DIR as argument or uses current directory
GPU_DIR="${1:-$(pwd)}"
cd "$GPU_DIR" || exit 1

DEVICE="${DEVICE:-CUDA}"       # Default device type
NUMWI="${NUMWI:-256}"           # Default work items
MAKE_J="${MAKE_J:-$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)}"


###autodock_gpu_<N>wi refers to the number of OpenCL/CUDA work-items per ligand 
#(like thread-group size):
#64wi = 64 threads per ligand
#128wi = 128 threads per ligand (more parallelism)

upper_device="$(echo "$DEVICE" | tr '[:lower:]' '[:upper:]')"
device_mode="CUDA"
if [[ "$upper_device" == "OPENCL" || "$upper_device" == "OCL" || "$upper_device" == "OCLGPU" || "$upper_device" == "AMD" ]]; then
  device_mode="OCLGPU"
fi
echo "[INFO] Using device mode: $device_mode"
# Standard CUDA binary naming
CUDA_BINARY="${GPU_DIR}/bin/autodock_gpu_cuda_${NUMWI,,}wi"
OCL_BINARY="${GPU_DIR}/bin/autodock_gpu_ocl_${NUMWI,,}wi" ##im not gonna use it but copilot suggested it, so lets keep it for now


if [[ "${DEVICE^^}" == "CPU" ]]; then
  echo "[INFO] CPU mode → skipping AutoDock-GPU build; proceeding to AutoGrid check/build…"
else
  ### Step 1: Compile AutoDock-GPU
  if [[ "$device_mode" == "CUDA" ]]; then
    if [ ! -f "$CUDA_BINARY" ]; then
        echo "[INFO] AutoDock-GPU binary not found. Compiling for $DEVICE..."

        if ! command -v nvcc &> /dev/null; then
            echo "[ERROR] CUDA toolkit not found. Please install it and try again."
            exit 1
        fi

        export GPU_INCLUDE_PATH="/usr/local/cuda/include"
        export GPU_LIBRARY_PATH="/usr/local/cuda/lib64"

        # Auto-detect TARGETS (keeps manual override if you set TARGETS/TARGET_ARCH)
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
        mv "${GPU_DIR}/bin/autodock_gpu_${NUMWI,,}wi" "${GPU_DIR}/bin/autodock_gpu_cuda_${NUMWI,,}wi"

        if [ -f "$CUDA_BINARY" ]; then
            echo "[INFO] AutoDock-GPU compilation successful."
        else
            echo "[ERROR] AutoDock-GPU compilation failed."
            exit 2
        fi
    else
        echo "[INFO] AutoDock-GPU binary already compiled."
    fi
  elif [[ "$device_mode" == "OCLGPU" ]]; then
    # OpenCL/AMD path
    echo "[INFO] Building AutoDock-GPU for OpenCL…"

    # Try to auto-fill headers/libs for OpenCL
    if [ -d /opt/rocm ]; then
      # ROCm path
      export GPU_INCLUDE_PATH="${GPU_INCLUDE_PATH:-/opt/rocm/include}"
      # libOpenCL.so often under /opt/rocm/lib or /opt/rocm/lib64
      if [ -f /opt/rocm/lib/libOpenCL.so ]; then
        export GPU_LIBRARY_PATH="${GPU_LIBRARY_PATH:-/opt/rocm/lib}"
      else
        export GPU_LIBRARY_PATH="${GPU_LIBRARY_PATH:-/opt/rocm/lib64}"
      fi
    else
      # Debian/Ubuntu ocl-icd-opencl-dev
      export GPU_INCLUDE_PATH="${GPU_INCLUDE_PATH:-/usr/include}"
      # Try to locate libOpenCL.so
      if ldconfig -p 2>/dev/null | grep -q 'libOpenCL\.so'; then
        OCL_LIB_DIR="$(ldconfig -p | awk '/libOpenCL\.so/{print $NF}' | head -n1 | xargs dirname)"
        export GPU_LIBRARY_PATH="${GPU_LIBRARY_PATH:-$OCL_LIB_DIR}"
      else
        export GPU_LIBRARY_PATH="${GPU_LIBRARY_PATH:-/usr/lib/x86_64-linux-gnu}"
      fi
    fi

    # Quick sanity for OpenCL toolchain
    if ! command -v clinfo >/dev/null 2>&1 && \
      [ ! -f "${GPU_INCLUDE_PATH}/CL/cl.h" ] && \
      ! ldconfig -p 2>/dev/null | grep -q 'libOpenCL\.so'
    then
      cat <<EOF
[ERROR] OpenCL headers/runtime not detected.
Install OpenCL ICD & headers. On Debian/Ubuntu:
  sudo apt-get update && sudo apt-get install -y ocl-icd-opencl-dev clinfo
For AMD ROCm, ensure /opt/rocm/* contains OpenCL headers/libs.
EOF
      exit 1
    fi

    # Already-built?
    OCL_PRESENT=""
    for name in autodock_gpu_ocl autodock_gpu_opencl autodock_gpu_ocl_${NUMWI}wi autodock_gpu_${NUMWI}wi_ocl autodock_gpu; do
      if [[ -x "${GPU_DIR}/bin/$name" ]]; then OCL_PRESENT="${GPU_DIR}/bin/$name"; break; fi
    done

    if [[ -n "$OCL_PRESENT" ]]; then
      echo "[INFO] OpenCL binary already present: $OCL_PRESENT"
    else
      mkdir -p ./bin
      echo "[INFO] Running make DEVICE=OCLGPU NUMWI=${NUMWI}…"
      make DEVICE=OCLGPU NUMWI="${NUMWI}" CC=gcc CXX=g++ -j"${MAKE_J}"
      export GPU_INCLUDE_PATH=/usr/include
      export GPU_LIBRARY_PATH=/lib/x86_64-linux-gnu
      # Quiet the CPU env checks (Makefile.OpenCL prints those)
      export CPU_INCLUDE_PATH=/usr/include
      export CPU_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
      mv "${GPU_DIR}/bin/autodock_gpu_${NUMWI,,}wi" "${GPU_DIR}/bin/autodock_gpu_ocl_${NUMWI,,}wi"


      OUT="${GPU_DIR}/bin/autodock_gpu_ocl_${NUMWI,,}wi"
      if [[ -x "$OUT" ]]; then
        echo "[INFO] OpenCL binary ready: $OUT"
      else
        echo "[ERROR] Build finished but $OUT not found"
        exit 2
      fi
    fi
  fi
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
    AUTOTOOLS_OK=true
    ACLOCAL_BIN=""
    CSH_OK=true
    if command -v aclocal-1.16 >/dev/null 2>&1; then
        ACLOCAL_BIN="$(command -v aclocal-1.16)"
    elif command -v aclocal >/dev/null 2>&1; then
        ACLOCAL_BIN="$(command -v aclocal)"
    else
        AUTOTOOLS_OK=false
        echo "[WARN] 'aclocal' not found; skipping autoreconf regeneration and suppressing automake hooks."
    fi

    if [ -n "$ACLOCAL_BIN" ]; then
        export ACLOCAL="$ACLOCAL_BIN"
    fi

    if ! command -v csh >/dev/null 2>&1; then
        CSH_OK=false
        echo "[WARN] 'csh' not found; will pre-generate default_parameters.h without it."
    fi

    if [ ! -f configure ]; then
        if [ -f configure.ac ] && [ "$AUTOTOOLS_OK" = true ]; then
            echo "[INFO] Running autoreconf to generate configure script..."
            autoreconf -i
        elif [ -f configure~ ]; then
            echo "[WARN] configure missing; restoring from configure~ backup."
            cp configure~ configure
        elif [ "$AUTOTOOLS_OK" = true ] && [ -f configure.in ]; then
            echo "[INFO] Found legacy configure.in; running autoreconf..."
            autoreconf -i
        else
            echo "[ERROR] configure script missing and no autoreconf sources available."
            echo "        Please re-fetch the AutoGrid sources or restore configure manually."
            exit 4
        fi
    fi

    if [ ! -x configure ] && [ -f configure ]; then
        chmod +x configure
    fi

    echo "[INFO] Running ./configure..."
    ./configure
    echo "[INFO] Cleaning previous build (if any)..."
    if [ "$AUTOTOOLS_OK" = true ]; then
        make clean || true
    else
        make clean ACLOCAL=: AUTOCONF=: AUTOMAKE=: AUTOHEADER=: || true
    fi

    if [ "$CSH_OK" = false ]; then
        python3 - "$AUTOGRID_DIR" <<'PY'
import sys
from pathlib import Path

autogrid_dir = Path(sys.argv[1])
out_file = autogrid_dir / "default_parameters.h"
srcs = [
    autogrid_dir / "ad4_shared" / "AD4_parameters.dat",
    autogrid_dir / "ad4_shared" / "AD4.1_bound.dat",
]

def load_lines(path: Path) -> list[str]:
    lines: list[str] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            if raw.lstrip().startswith("#"):
                continue
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            lines.append(line)
    return lines

blocks = [
    ("const char *param_string_4_0[MAX_LINES] = {", load_lines(srcs[0])),
    ("const char *param_string_4_1[MAX_LINES] = {", load_lines(srcs[1])),
]

with open(out_file, "w", encoding="utf-8") as out:
    for header, entries in blocks:
        out.write(f"{header}\n")
        for entry in entries:
            out.write(f"\"{entry}\\n\", \n")
        out.write(" };\n")
    out.write("// EOF\n")
PY
        # ensure the header timestamp is newer than its dependencies so make skips the csh rule
        touch "$AUTOGRID_DIR/default_parameters.h"
        echo "[INFO] default_parameters.h generated via Python fallback."
    fi
    echo "[INFO] Running make..."
    if [ "$AUTOTOOLS_OK" = true ]; then
        make -j$(nproc)
    else
        make -j$(nproc) ACLOCAL=: AUTOCONF=: AUTOMAKE=: AUTOHEADER=:
    fi

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
