# Installation and screening preparation

## Choose a machine and backend

Ultidock requires Python 3.10 or later and a POSIX environment. Linux is the
documented native build path. CPU mode runs AutoDock Vina; CUDA and OpenCL modes
run AutoDock-GPU. AutoGrid is used for grid/site preparation in both paths. A
Python package installation does not install native engines or GPU drivers.

| Planned run | Hardware and software to prepare | Expected effect |
| --- | --- | --- |
| Small pilot or CPU screen | Linux host, Vina, AutoGrid build tools, enough CPU cores and RAM for concurrent jobs | Works without a GPU. More workers can process independent ligand/site jobs, while each Vina job uses `--vina-cpu` threads. |
| NVIDIA screen | Supported NVIDIA GPU, working driver and compatible CUDA toolkit, plus the native build tools | AutoDock-GPU searches poses; CPU preparation, site finding, I/O and analysis still take time. |
| OpenCL screen | Supported AMD/Intel/OpenCL GPU, vendor or Mesa device runtime, OpenCL headers and linker library | AutoDock-GPU uses the OpenCL backend when a real GPU device is visible; an ICD loader alone is insufficient. |
| Multiple sites or larger boxes | Any backend, with more RAM and free disk for grids and poses | Work grows with ligand–site pairs; more/larger grids also increase setup time and storage. |

The README recommends **16 GB RAM and 20 GB free storage for large batches** as
planning starting points, not fixed minimums. Actual requirements depend on the
receptor, box size, site count, ligand library and retained artifacts. Keep the
working directories on a disk with room for grids, raw poses, logs and the SQLite
database. See [Measure and tune a pilot](../tutorials/virtual-screening.md#5-measure-and-tune-a-pilot)
before reserving a large machine or cluster allocation.

## Install Python and native packages

From a source checkout, the documented **Ubuntu 26.04** package set covers the
CPU workflow, native builds, raw PDB conversion and the optional pocket tools:

```bash
sudo apt update
sudo apt install -y \
  autoconf automake autodock-vina build-essential clinfo cmake csh curl \
  g++-12 gcc-12 gfortran git libnetcdf-dev libtool libx11-dev \
  m4 make ocl-icd-opencl-dev openbabel openjdk-21-jre-headless \
  perl pkg-config python3 python3-pip python3-venv tar unzip wget
```

`autodock-vina` provides `vina` and `vina_split`; Open Babel supports receptor
conversion and donor-hydrogen recovery. The listed compilers and build tools
cover the included AutoGrid/AutoDock-GPU sources; `libnetcdf-dev` is used when
building fpocket, and Java 21 is the documented choice for P2Rank 2.5. Other
Linux distributions need equivalent packages. `requirements.txt` contains only
Python packages.

Install Ultidock and its Python dependencies from the repository root:

```bash
git clone --branch v1.1.2 https://github.com/taka78/ultidock.git
cd ultidock
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
ultidock --help
ultidock doctor
```

Keep the venv activated for later commands. An editable development install is
`python -m pip install -e ".[dev]"`. With Ubuntu system Python packages, you can
instead use `/usr/bin/python3 -m cli.ultidock` and
`/usr/bin/python3 -m cli.molguard` from the checkout; the packages must satisfy
`pyproject.toml`. Do not mix that interpreter with packages from a different venv.

For raw PDB receptor conversion, the package command above installs Open Babel.
Meeko is another receptor converter and can prepare ligands from suitable source
structures. If you choose it, install it in the active Python environment with
`python -m pip install meeko` and check its preparation commands. Prepared
receptor PDBQT files can be validated without converting a raw PDB. See
[Preparing a receptor](../user-guide/preparing-receptor.md) and
[Preparing ligands](../user-guide/preparing-ligands.md).

## Prepare a GPU runtime when needed

For CUDA, install a driver and CUDA toolkit compatible with the GPU and native
AutoDock-GPU build; verify `nvidia-smi -L` and `nvcc --version`. For OpenCL,
install a device runtime as well as the headers/linker library and verify
`clinfo -l` lists a GPU device. The Ubuntu package command supplies OpenCL
development files, not a vendor runtime. For AMD Mesa/Rusticl, the repository
[GPU runtime notes](https://github.com/taka78/ultidock/blob/v1.1.2/README.md#gpu-runtimes)
give the package and device checks.

```bash
ultidock setup --mode cpu --skip-wget   # use Vina; build/reuse AutoGrid
# or, after verifying the GPU runtime:
ultidock setup --mode gpu --skip-wget   # AutoDock-GPU and AutoGrid
ultidock doctor
```

The first setup may build native tools. `--mode gpu` stops when detection fails;
`--mode auto` permits CPU fallback. Select an explicit mode for a measured screen
so the engine does not change between runs. Device detection alone does not prove
that a complete docking calculation succeeds; run the bundled example next.

## Optional local pocket predictors

CaV-EMPS is built into the default `cavity` workflow. Install fpocket or P2Rank
only when selecting those methods. From the checkout:

```bash
bash scripts/install_pocket_tools.sh fpocket
bash scripts/install_pocket_tools.sh p2rank
ultidock doctor
```

The installers target fpocket 4.2.3 and P2Rank 2.5. P2Rank needs Java 17–23;
Java 21 is the documented choice. A regular installation can also install the
selected local predictor when first invoked. On a cluster, provision it before
submitting many jobs. See the [binding-site guides](../user-guide/binding-site-discovery/index.md).

## Verify a real run and locate the workspace

```bash
ultidock example run d2-antipsychotics --dry-run
ultidock example run d2-antipsychotics --mode cpu
```

The preview does not stage or dock; the real example stages one receptor and
three ligands in a timestamped workspace. Review the path printed by the runner
and check for poses and results there. For a GPU validation, repeat with
`--mode gpu` after its setup succeeds. Then follow
[Virtual screening](../tutorials/virtual-screening.md) for your own library.

For a regular installation, writable workflow files live under
`$XDG_DATA_HOME/ultidock` (default `~/.local/share/ultidock`); set
`ULTIDOCK_HOME=/path/to/workspace` to choose another location. Editable and
direct checkout runs use the repository layout by default. `ultidock doctor`
prints the active workspace. Keep separate homes for concurrent runs because
setup writes `docking/config.py` there.
