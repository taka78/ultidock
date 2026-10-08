# System requirements and GPU runtimes

Ultidock's native build path is Linux. Ubuntu 26.04 is the documented package
example; other Linux distributions need equivalent tools and runtimes. Windows
and macOS users should use a Linux container or virtual machine with access to
the needed devices. A Bash shell and coreutils should be on `PATH`.

## Hardware planning

| Component | Starting point |
| --- | --- |
| CPU | x86-64 with AVX for preprocessing and optional Vina docking |
| GPU | Supported NVIDIA CUDA or OpenCL device for AutoDock-GPU; CPU mode needs no GPU |
| RAM | 16 GB recommended for large ligand batches |
| Storage | At least 20 GB free as a planning starting point for archives, grids and outputs |

These are planning figures, not guaranteed minima.
Receptor size, box dimensions, site count, ligand count, worker count and
retained raw artifacts determine actual memory, VRAM, time and disk use.
Measure a representative [screening pilot](../tutorials/virtual-screening.md)
on the machine that will run the full library.

## Ubuntu native packages

Install the tools for AutoGrid, AutoDock-GPU builds, CPU Vina, fpocket, P2Rank
and raw PDB receptor conversion:

```bash
sudo apt update
sudo apt install -y \
  autoconf automake autodock-vina build-essential clinfo cmake csh curl \
  g++-12 gcc-12 gfortran git libnetcdf-dev libtool libx11-dev \
  m4 make ocl-icd-opencl-dev openbabel openjdk-21-jre-headless \
  perl pkg-config python3 python3-pip python3-venv tar unzip wget
```

`autodock-vina` supplies `vina` and `vina_split`; Open Babel converts raw
receptors and supports donor-hydrogen recovery. AutoGrid's build needs
Autotools, `m4`, Perl and `csh`. The current AutoDock-GPU build script calls
`gcc-12` and `g++-12`. fpocket needs a C/C++ toolchain and NetCDF headers;
P2Rank needs `curl`, `tar` and Java 17–23 (Java 21 is the documented choice).
`requirements.txt` installs Python packages only. A Python environment does
not supply a driver or vendor compute runtime.

## NVIDIA CUDA

Install a driver and CUDA toolkit compatible with the GPU and the AutoDock-GPU
build. Check from the same terminal, container or compute node that runs
Ultidock:

```bash
nvidia-smi -L
nvcc --version
```

The first command must list a device and the second must find the compiler for
a fresh CUDA build. Consult [NVIDIA's CUDA downloads](https://developer.nvidia.com/cuda-downloads)
for a compatible installation. A working desktop display by itself does not
establish CUDA compute access. Containers and VMs need device and runtime
access too.

## AMD or Intel OpenCL

Install a supported vendor or Mesa device runtime, the OpenCL ICD development
headers/linker library and `clinfo`. Check that `clinfo -l` lists a **GPU
device**, not merely a platform. `Number of platforms 0` means no OpenCL
platform is visible to that process. A visible device still needs a successful
AutoDock-GPU build and test run.

On Ubuntu with AMD Mesa/Rusticl, the documented package and visibility check is:

```bash
sudo apt install mesa-opencl-icd ocl-icd-opencl-dev clinfo
export RUSTICL_ENABLE=radeonsi
clinfo -l
```

`mesa-opencl-icd` supplies the Mesa runtime;
`ocl-icd-opencl-dev` supplies `CL/opencl.h`, `CL/cl.h` and
`libOpenCL.so` for compilation. A working `clinfo` does not prove the
development files are installed. The manual `RUSTICL_ENABLE` export is useful
for checking the device. See [Mesa's environment variable documentation](https://docs.mesa3d.org/envvars.html#envvar-RUSTICL_ENABLE).

If Rusticl reports no GPU and `RUSTICL_ENABLE` is unset, Ultidock setup retries
with `radeonsi`. It selects OpenCL only if the retry exposes a GPU, passes the
setting to subprocesses, and stores it in `docking/config.py` for later
`--skip-setup` runs. It respects an existing user setting and does not edit a
shell profile. This applies to explicit `--mode opencl` setup too.

## Backend selection

The default `--mode auto` prefers a detected NVIDIA CUDA GPU, then a visible
OpenCL GPU, and uses CPU Vina when none is detected. On CUDA, the runner
distributes jobs across detected NVIDIA GPU IDs, with two slots per device by
default. Check the logged GPU IDs and memory use when scaling up.
`--mode gpu` requires detection and stops when no GPU is found.
`--mode cpu` selects Vina without a GPU.
`--mode cuda` and `--mode opencl` select those backends explicitly.
AutoGrid is built or reused even in CPU mode. Use an explicit mode for an
experiment so its engine stays fixed, then validate it with a complete
[first docking run](first-docking-run.md).

See [Installation](installation.md) for Python setup and
[Troubleshooting](../reference/troubleshooting.md) for build failures.
