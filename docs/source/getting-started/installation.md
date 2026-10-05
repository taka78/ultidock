# Installation

## Supported starting point

The workflow expects a POSIX environment; Linux is the documented native build
path. Python 3.10 or later is required. The Python package provides `ultidock` and
`molguard`; native engines and GPU runtimes are separate dependencies.

From a checkout:

```bash
git clone --branch v1.1.2 https://github.com/taka78/ultidock.git
cd ultidock
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
ultidock --help
ultidock doctor
```

Keep the environment activated in subsequent commands. As an alternative, a
checkout can run with `/usr/bin/python3 -m cli.ultidock` and
`/usr/bin/python3 -m cli.molguard` when all dependencies in `pyproject.toml` are
installed for system Python. Do not mix that interpreter with a different venv's packages.

## Native tools

On Ubuntu, a CPU-oriented starting point is:

```bash
sudo apt update
sudo apt install autodock-vina openbabel build-essential autoconf automake \
  libtool csh git wget curl unzip
```

Vina performs CPU docking. AutoGrid is also used by the main Ultidock grid/site
pipeline. Setup reuses available tools or attempts to build the bundled native
sources. Consult the repository's
[complete setup guide](https://github.com/taka78/ultidock/blob/v1.1.2/SETUP.md)
for the full compiler/build dependencies, including the GCC toolchain expected
by AutoDock-GPU and optional pocket tools.

For Meeko-based receptor conversion, install its receptor dependencies in your
chosen Python environment; an available `mk_prepare_receptor.py` is detected by
MolGuard. Open Babel is an alternative and is required for the PDBQT donor-hydrogen
recovery path. See [Preparing a receptor](../user-guide/preparing-receptor.md).

## GPU and pocket methods

CUDA requires a compatible NVIDIA driver/runtime; OpenCL requires an ICD and a
supported device runtime. `--mode gpu` fails if no suitable device is found;
`--mode auto` explicitly permits CPU fallback. A Python package install cannot
supply these vendor drivers.

Install optional pocket tools from the checkout:

```bash
bash scripts/install_pocket_tools.sh fpocket
bash scripts/install_pocket_tools.sh p2rank
ultidock doctor
```

P2Rank's bundled installer targets 2.5 and a Java 17–23 runtime; Java 21 is the
documented choice. fpocket needs its C/build dependencies. See the method guides
before using either tool on a cluster.

Finish by previewing the [Quick Start](quick-start.md):

```bash
ultidock example run quickstart --dry-run
```
