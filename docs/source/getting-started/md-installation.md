# Install the molecular dynamics tools

The `gmx-dev` branch adds `ultidock md`, which takes scored docking poses into
GROMACS molecular dynamics (MD). Install these additions in the same environment
as Ultidock. Docking setup is covered by [Install Ultidock](installation.md).

## Required tools

| Tool | Purpose | Availability check |
| --- | --- | --- |
| GROMACS 2021 or newer | Build the system and run minimization, equilibration and production | `gmx --version` |
| AmberTools | GAFF2 atom types and AM1-BCC ligand charges | `antechamber`, `parmchk2`, `tleap`, `sqm` on `PATH` |
| ACPYPE | Convert Amber ligand parameters to GROMACS topology | `acpype` on `PATH` |
| Open Babel | Convert the chemically complete posed ligand to MOL2 | `obabel` and Python bindings for ACPYPE |
| RDKit, NumPy, SciPy | Ligand chemistry, atom mapping and system assembly | Modules in Ultidock's Python environment |

Use a **thread-MPI `gmx` executable**. The current runner starts GROMACS
directly and sets its thread budget; it does not launch external MPI jobs.
GROMACS 2021+ is needed for the workflow's C-rescale pressure coupling.

## Native packages and your existing Python environment

On Ubuntu/Debian, install GROMACS and Open Babel, then check the GROMACS version:

```bash
sudo apt update
sudo apt install gromacs openbabel
gmx --version
```

From the Ultidock checkout root, activate your Python environment. If needed,
create one first:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[md]" openbabel-wheel
```

The `md` extra installs RDKit and ACPYPE in addition to the ordinary Ultidock
dependencies. Keep ACPYPE and its Open Babel bindings in that environment;
see the [ACPYPE installation guide](https://github.com/alanwilter/acpype/blob/main/README.md).
Available ACPYPE wheels and their bundled native tools depend on the platform.
Ultidock checks the four AmberTools commands separately, so verify their paths
with `ultidock doctor` rather than assuming a pip installation supplied them.

## Standalone AmberTools

If AmberTools is already installed, source its `amber.sh`. Otherwise download
the source distribution from [Amber's download page](https://ambermd.org/GetAmber.php#ambertools)
and follow its reference manual and [Ubuntu build instructions](https://ambermd.org/InstUbuntu.php).
The repository's documented build also needs `cmake` and `python3-tk`.
Install the Python build dependencies in your active environment:

```bash
python -m pip install numpy scipy matplotlib cython setuptools
python -c 'import sys; print(sys.executable)'
```

Before running the source tree's `build/run_cmake`, set these CMake options
using the interpreter path printed above and your intended installation prefix:

```text
-DDOWNLOAD_MINICONDA=FALSE
-DPYTHON_EXECUTABLE=/absolute/path/to/ultidock/.venv/bin/python
-DCMAKE_INSTALL_PREFIX=/absolute/path/to/ultidock/md-simulation/tools/ambertools
```

These [Amber CMake options](https://ambermd.org/pmwiki/pmwiki.php/Main/CMake-Common-Options)
use your existing Python environment and disable the build's Miniconda download.
Inside the extracted AmberTools source tree's `build` directory:

```bash
./run_cmake
cmake --build . --parallel 4
cmake --install .
```

Return to the Ultidock root and load the installation:

```bash
source md-simulation/tools/ambertools/amber.sh
ultidock doctor
```

Use your actual prefix if different. In a new terminal, activate the Python
environment and source `amber.sh` again before running Ultidock.

## Optional Conda environment

The root README also documents a complete environment under `md-simulation/tools/`:

```bash
conda create --prefix ./md-simulation/tools/env \
  --override-channels --channel conda-forge --strict-channel-priority \
  python=3.12 pip "gromacs=*=nompi_h*" ambertools acpype openbabel rdkit
conda activate ./md-simulation/tools/env
python -m pip install -e .
ultidock doctor
```

The GROMACS selector requests the CPU, non-MPI build that supplies `gmx`.
Check the [conda-forge build definitions](https://github.com/conda-forge/gromacs-feedstock/blob/main/recipe/meta.yaml)
if that selector is unavailable for your platform. Activate the same environment
in later terminals using its full path.

## Custom GROMACS and executable paths

For a source-built or GPU-enabled GROMACS, follow the
[GROMACS installation guide](https://manual.gromacs.org/current/install-guide/index.html)
and source the installed `bin/GMXRC`. GPU acceleration depends on the GROMACS
build and supported runtime; installing AutoDock-GPU alone does not enable it.

```bash
ultidock doctor --gmx /path/to/gmx --acpype /path/to/acpype --obabel /path/to/obabel
```

Use those same overrides with `ultidock md run`. Integrated docking uses the
corresponding `--md-gmx`, `--md-acpype` and `--md-obabel` options. AmberTools
commands must remain on `PATH`.

Doctor should report `MD dependencies: ready`. It checks availability, while
the [docking-to-MD workflow](../user-guide/molecular-dynamics/docking-to-md.md)
checks the protocol, molecular inputs and native simulation stages. See
[running and results](../user-guide/molecular-dynamics/running-and-results.md)
for thread settings, GPU behavior and performance planning.
