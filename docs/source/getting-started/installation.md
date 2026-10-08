# Install Ultidock

Install the [native packages and any GPU runtime](requirements.md) first.
Ultidock requires Python 3.10 or later. The package installs the
`ultidock` workflow CLI, the `molguard` validation CLI, and Python
dependencies from `pyproject.toml`: Click, NumPy, SciPy, psutil, pandas and
Matplotlib (plus `tomli` on Python before 3.11). Python packages alone do
not install docking engines, compilers, Java or GPU drivers.

## Recommended: virtual environment

From a checkout:

```bash
git clone --branch gmx-dev https://github.com/taka78/ultidock.git
cd ultidock
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
ultidock --help
```

Activate `.venv` again in each new terminal. There is no `make install`
target in the repository root. For development, use
`python -m pip install -e ".[dev]"`. If upgrading an environment with the
old separate `molguard` distribution, remove it first with
`python -m pip uninstall molguard`, then install Ultidock. The
`molguard` command and Python imports remain part of Ultidock.

## Alternative: Ubuntu system Python

On Ubuntu 26.04, install sufficiently recent distro Python packages without
a venv or pip installation of this checkout:

```bash
sudo apt install python3 python3-click python3-numpy python3-scipy \
  python3-psutil python3-pandas python3-matplotlib
/usr/bin/python3 -m cli.ultidock --help
/usr/bin/python3 -m cli.ultidock example run quickstart
/usr/bin/python3 -m cli.molguard --help
```

Run these commands from the repository root. `/usr/bin/python3` selects
system Python even if a venv is active. For this route, substitute
`/usr/bin/python3 -m cli.ultidock` for `ultidock` and
`/usr/bin/python3 -m cli.molguard` for `molguard` throughout the guides.
On other distributions, the installed versions must meet `pyproject.toml`.
Do not mix this interpreter with packages from another venv.

## Native tools and receptor conversion

`vina` and `vina_split` must be on `PATH` or available in the checkout's
bundled Linux Vina directory. The Ubuntu command in
[System requirements](requirements.md) installs `autodock-vina`.
An available compatible `autogrid4` is reused; otherwise setup builds it locally.
For CaV-EMPS, a binary advertising fewer than 20 maps is rebuilt from the bundled
sources. `ultidock doctor` reports this limitation.
Regular installs link available native programs into their writable workspace.
The wheel includes build **sources**, not machine-specific executables.

For raw PDB or MOL2 receptor conversion, use an installed Meeko receptor
command or Open Babel. The Ubuntu native-package command installs Open Babel;
an alternative in an activated venv is:

```bash
python -m pip install meeko
```

Prepared PDBQT receptors can be canonicalized without a converter. This branch
does not add missing donor hydrogens to an existing PDBQT; prepare it again from
a chemically appropriate source if needed.
Meeko can also prepare ligands from suitable source structures. Choose
chemistry and conformers before conversion. See
[Preparing a receptor](../user-guide/preparing-receptor.md) and
[Preparing ligands](../user-guide/preparing-ligands.md).

## Managed workspace

A regular installation can run outside the source checkout. On first workflow
use, Ultidock materializes scripts, examples, benchmark code and native build
sources under `$XDG_DATA_HOME/ultidock`, defaulting to
`~/.local/share/ultidock`. Configuration, compiled tools, predictor
downloads and run results live there; installed Python files remain untouched.
Set `ULTIDOCK_HOME=/path/to/workspace` to choose another location.
Package upgrades refresh managed source files while retaining configuration
and results. Editable and direct checkout runs use the repository layout by
default, but also honor `ULTIDOCK_HOME`. `ultidock doctor` prints the
active workspace.

## Optional fpocket and P2Rank

CaV-EMPS is included with Ultidock. Install only the alternative pocket tools
you plan to use:

```bash
bash scripts/install_pocket_tools.sh fpocket
bash scripts/install_pocket_tools.sh p2rank
ultidock doctor
```

Run the scripts from the checkout. With no argument or `all`, the installer
installs both. It targets fpocket 4.2.3 and P2Rank 2.5 under the ignored
`external/` directory. fpocket needs its C toolchain and NetCDF headers;
the installer also corrects a GCC 15 pointer-type error in the pinned source.
P2Rank needs Java 17–23, with Java 21 the documented choice.
Missing optional tools do not block setup; selecting a pocket method installs
only that tool if absent. `--tool PATH` can choose a preinstalled executable.
Provision tools once before launching a cluster array.

For Java incompatibility, see [P2Rank](../user-guide/binding-site-discovery/p2rank.md)
and [Troubleshooting](../reference/troubleshooting.md).

## Verify the installation

```bash
ultidock doctor
ultidock example run quickstart --dry-run
ultidock example run sert-escitalopram --dry-run
```

`doctor` distinguishes a tool whose source exists but has not been compiled
from a tool that cannot be found. Dry runs check dispatch and inputs without
proving the native engines work. Follow
[Dock your own molecules](../user-guide/start-docking.md) and complete a real
[first docking run](first-docking-run.md) before a large screen.

## Optional molecular dynamics

Docking alone uses the base installation above. To continue into MD, install
GROMACS, AmberTools, ACPYPE, RDKit and Open Babel following
[MD installation](md-installation.md). Then complete the reviewed protein,
ligand and system inputs in the [MD guide](../user-guide/molecular-dynamics/index.md).
`--md-config` connects the resulting protocol to a docking run without requiring
separate manual pose preparation.
