# Ultidock (gmx-dev branch)

Ultidock is a high-throughput molecular docking workflow that automates ligand
staging, grid preparation, AutoDock-GPU execution, and post-processing. The
`gmx-dev` branch focuses on reproducible automation today and prepares the
groundwork for future GROMACS-based molecular dynamics integration.

This document explains **how to run the pipeline step by step**, details the
major components, and highlights the features that make Ultidock different from
traditional docking scripts.

---

## Table of Contents

1. [Requirements](#requirements)
2. [Repository Layout](#repository-layout)
3. [Quick Start: End-to-End Run](#quick-start-end-to-end-run)
4. [Command Reference](#command-reference)
5. [Configuration Reference](#configuration-reference)
6. [Pipeline Segments & What Makes Ultidock Different](#pipeline-segments--what-makes-ultidock-different)
7. [Spotlight: Grid Boxing & Cavity Finder Algorithm](#spotlight-grid-boxing--cavity-finder-algorithm)
8. [Working with the Example Pipelines](#working-with-the-example-pipelines)
9. [Troubleshooting](#troubleshooting)
10. [Citation & License](#citation--license)

---

## Requirements

Ultidock targets modern Linux systems. Windows and macOS users should rely on a
Linux container or VM.

### Hardware

| Component | Requirement |
|-----------|-------------|
| CPU       | x86-64 with AVX (for preprocessing and optional CPU docking) |
| GPU       | NVIDIA GPU with CUDA capability **7.0 or newer** (Ampere, Ada, Hopper, or RTX 40/50). CPU-only mode is supported but slower. |
| RAM       | ≥ 16 GB recommended for large ligand batches |
| Storage   | ≥ 20 GB free space for ligand archives, grids, and outputs |

### Operating System

- Ubuntu 22.04+, Debian 12+, Fedora 39+, or a comparable modern Linux distro
- Bash shell and coreutils available on `$PATH`

### System Packages

Install the build toolchain and helper utilities once:

```bash
sudo apt update && sudo apt install -y \
  automake autoconf libtool m4 perl pkg-config\
  build-essential gcc g++ gfortran make cmake \
  unzip tar csh wget git \
  libstdc++-dev libx11-dev libncurses-dev \
  python3 python3-venv python3-pip
```

> **Tip:** Replace `apt` commands with the equivalent package manager commands
> for your distribution.

### GPU Runtimes

- **Latest available CUDA Toolkit for your hardware** is required for NVIDIA GPU execution. Install
  it from [NVIDIA's official downloads](https://developer.nvidia.com/cuda-downloads).
- Ultidock defaults to AutoDock-GPU. AutoGrid will also be compiled on first run.

### Python Environment

Ultidock requires Python **3.10+**. A virtual environment is recommended:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
```

Install all required packages and the `ultidock` CLI in one step:

```bash
pip install -r requirements.txt
pip install -e .
```

This installs numpy, scipy, psutil, matplotlib, pandas and the `molguard`
I/O hardening layer (exposed as the `ultidock` command).
`pandas` and `matplotlib` are used only by the post-run analysis stage.

---

## Repository Layout

```
ultidock/
├─ pyproject.toml            # Package manifest; defines the `ultidock` CLI entry-point
├─ requirements.txt          # Runtime dependencies (numpy, scipy, click, ...)
├─ requirements-dev.txt      # Development/test dependencies
├─ SETUP.md                  # Step-by-step new-user guide
├─ molguard/                 # I/O hardening + validation layer
│  ├─ io/
│  │  ├─ fixedfmt.py         # Fixed-width, locale-safe float formatter
│  │  └─ pdbqt.py            # PDBQT linter, normalizer, receptor canonicalizer
│  ├─ grids/
│  │  └─ check.py            # AutoGrid .fld / .map sanity checker
│  ├─ cli.py                 # `ultidock` command-line interface
│  └─ tests/                 # Unit and regression tests (pytest)
├─ docking/
│  ├─ run.py                 # Main entry point for the entire pipeline
│  ├─ setup.py               # Idempotent environment + dependency setup
│  ├─ dock_v02.py            # AutoDock-GPU / AutoGrid orchestration
│  ├─ analyse_docking_results.py
│  ├─ extract.py             # Front-end for AutoDock Vina's vina_split utility
│  ├─ clean.py               # Resets compiled binaries and outputs
│  ├─ ligands.wget           # Example ligand download manifest
│  ├─ MACRO_MOL_DIR/         # Receptors, grids, and generated sites
│  ├─ LIGANDS_DIR/           # Archived ligands and split PDBQT files
│  ├─ DOCKING_DIR/           # AutoDock-GPU/Vina output poses
│  ├─ AUTODOCK_GPU_DIR/      # Compiled AutoDock-GPU + AutoGrid binaries
│  ├─ VINA_DIR/              # AutoDock Vina binaries
│  ├─ ANALYSIS_DIR/          # Intermediate scoring/aggregation artifacts
│  └─ RESULTS_DIR/           # Final CSV/JSON summaries
├─ examples/                 # Self-contained example runners
└─ data-analyses/, results/  # Optional downstream notebooks & exports
```

The workflow assumes you copy or generate receptor `.pdbqt` files inside
`docking/MACRO_MOL_DIR/` and provide a `.wget` manifest (or existing ligand
archives) inside `docking/LIGANDS_DIR/`.

---

## Quick Start: End-to-End Run

Follow this checklist whenever you want to run Ultidock from a clean workspace.

1. **Clone the repository** (or update your local copy):
   ```bash
   git clone https://github.com/taka78/ultidock.git
   cd ultidock
   ```

2. **Activate your Python environment** and install the requirements (see
   [Requirements](#requirements)).

3. **Reset the docking workspace** to avoid stale binaries and outputs:
   ```bash
   python3 docking/clean.py -y --all
   ```

4. **Stage inputs:**
   - Copy your receptor(s) to `docking/MACRO_MOL_DIR/`. Each receptor can live in
     its own subdirectory if you plan to run multi-site docking.
   - Provide ligands via one of the following:
     - Populate `docking/ligands.wget` with direct links to `.pdbqt.gz` archives
       (one per line). Ultidock will download, verify, and extract them.
     - Manually place `.pdbqt` or `.pdbqt.gz` files in `docking/LIGANDS_DIR/`.
     - Pass `--skip-wget` when running `setup.py`/`run.py` to skip downloads and
       rely entirely on pre-populated ligand files.

5. **Run the setup + docking pipeline:**
   ```bash
   python3 docking/run.py --mode gpu
   ```
   - Use `--mode cpu` to skip AutoDock-GPU compilation and rely on AutoGrid
     + Vina.
   - Add `--skip-wget` if your ligands are already staged and you want to avoid
     executing the download manifest.
   - Override directories as needed with `--LIGANDS_DIR`, `--MACRO_MOL_DIR`, etc.
     Absolute paths are recommended for scripted automation.

6. **Monitor progress:**
   - Setup output reports where AutoDock-GPU, AutoGrid, and Vina binaries are
     compiled or reused.
   - Docking output prints the number of ligands discovered, grid preparation
     steps, worker launches, and database insertions.

7. **Review results:**
   - Raw poses are written to `docking/DOCKING_DIR/`.
   - Per-receptor metadata (centers, grids) lives in `docking/MACRO_MOL_DIR/`.
   - An auto-managed SQLite database (`docking/RESULTS_DIR/ultidock_results.db`)
     is updated throughout the run for incremental result parsing and can be
     inspected or queried at any time.
   - If `pandas` is installed, aggregated CSV/JSON summaries will be produced in
     `docking/RESULTS_DIR/`.

8. **Optional post-run steps:**
   - Run `python3 docking/extract.py --help` to (re)split ligand archives via
     AutoDock Vina's `vina_split` utility or prepare filtered subsets for
     downstream MD.
   - Use the notebooks in `data-analyses/` for visualization or scoring audits.

Repeat steps 3–8 for each new batch to ensure deterministic runs.

---

## Command Reference

| Command | Purpose |
|---------|---------|
| `python3 docking/run.py [options]` | Primary entry point. Validates the environment, runs setup, downloads ligands, launches docking, and triggers analysis. |
| `python3 docking/setup.py [options]` | Runs the setup stage only (directory creation, AutoDock-GPU/AutoGrid/Vina checks). All CLI flags mirror `run.py`. |
| `python3 docking/dock_v02.py [options]` | Executes the docking stage against prepared ligands and receptors. Used internally by `run.py`. |
| `python3 docking/extract.py` | Wrapper around AutoDock Vina's `vina_split` for splitting ligand archives and optional filtering. |
| `python3 docking/clean.py -y --all` | Removes compiled binaries, cached grids, downloads, and generated configs. Use before starting a fresh run. |
| `ultidock pdbqt check <file>` | Lint a receptor or ligand PDBQT for AutoDock column-format issues (exponent notation, missing decimals, bad atom types). |
| `ultidock pdbqt normalize <file> -o <out>` | Rewrite all numeric columns in a ligand PDBQT through the fixed-width formatter. Torsion tree is left untouched. |
| `ultidock pdbqt canonicalize-receptor <file> -o <out>` | Sort, renumber, and reformat a receptor PDBQT deterministically. Returns a SHA-256 digest for reproducibility checks. |
| `ultidock grids check <maps.fld>` | Validate AutoGrid output: checks for all-zero maps, NaN/Inf energies, missing files, and atom-type mismatches. |
| `ultidock doctor` | Print tool locations and versions. Distinguishes between binaries not compiled yet (source present) and not found at all. |

### Key `run.py` Flags

| Flag | Description |
|------|-------------|
| `--mode {gpu,cpu}` | Select GPU (AutoDock-GPU) or CPU-only (Vina) execution mode. |
| `--skip-setup` | Assume setup has already been run and use the existing config. |
| `--LIGANDS_DIR PATH` | Override ligand staging directory. |
| `--MACRO_MOL_DIR PATH` | Override receptor directory. |
| `--AUTODOCK_GPU_DIR PATH` | Override AutoDock-GPU build/install directory. |
| `--VINA_DIR PATH` | Override AutoDock Vina install directory. |
| `--RESULTS_DIR PATH`, `--ANALYSIS_DIR PATH`, `--DOCKING_DIR PATH` | Customize other pipeline locations. |
| `--wget FILE` | Use a custom `.wget` manifest for ligand downloads. |
| `--skip-wget` | Skip executing `wget` commands even if a manifest is present. |

All flags are optional; defaults point to directories within `docking/`.

---

## Configuration Reference

Running `python3 docking/run.py` or `python3 docking/setup.py` writes a fully
resolved configuration to `docking/config.py`. The file records the exact
directories, binaries, and grid parameters that Ultidock will reuse on the next
invocation. Edit the file directly (or pass CLI overrides) to fine-tune a run.

### Directory Layout Variables

| Variable | Meaning |
|----------|---------|
| `LIGANDS_DIR` | Absolute path where ligand archives and split PDBQT files are staged. |
| `DOCKING_DIR` | Output directory for AutoDock-GPU / Vina poses and logs. |
| `ANALYSIS_DIR` | Workspace for intermediate scoring, per-ligand summaries, and temporary exports. |
| `VINA_DIR` | Location of the AutoDock Vina binaries used for ligand splitting or CPU docking. |
| `AUTODOCK_GPU_DIR` | Location of the AutoDock-GPU and AutoGrid toolchains compiled during setup. |
| `MACRO_MOL_DIR` | Root folder for receptor structures, generated grids, and per-site artifacts. |
| `RESULTS_DIR` | Destination for final CSV/JSON exports and the SQLite results database. |
| `DB_PATH` | Full path to the SQLite database (`ultidock_results.db`) that receives live docking updates. |

### Runtime Controls

| Variable | Description |
|----------|-------------|
| `GPU_TYPE` | Which accelerator build to prepare (`CPU`, `CUDA`, or `OCL`). In CPU mode only AutoGrid and Vina are compiled. |
| `NUMWI` | Number of AutoDock-GPU work items queued per ligand batch. Increase to better saturate large GPUs; reduce on memory-constrained devices. |
| `AUTO_GRID_BIN` | Resolved path to the `autogrid4` binary. Adjust if you provide a prebuilt AutoGrid installation. |
| `GRID_MODE` | Strategy for identifying grid centers: `ligand`, `residues`, `centers` (hotspot-driven default), or `blind` (whole-protein). |
| `GRID_SPACING` | Ångström spacing between grid points. Smaller values yield finer resolution at the cost of longer AutoGrid runtimes. |
| `GRID_MARGIN` | Extra Ångström padding applied to each hotspot-derived grid to ensure the box fully encloses the binding site. |
| `GRID_CAP` | Maximum Å-length per axis when running in blind mode to prevent runaway grid sizes. |
| `CENTERS_TSV` | Optional path to a precomputed `centers.tsv`. Leave as `None` to let Ultidock regenerate hotspot centers automatically. |
| `REF_LIGAND_PDB` | Reference ligand file used when `GRID_MODE="ligand"` to seed the search box from a co-crystal pose. |

### Grid Boxing & Cavity Finder Dials

These parameters feed the hotspot detection routine acknowledged in the
[Spotlight](#spotlight-grid-boxing--cavity-finder-algorithm) section.

| Variable | Description |
|----------|-------------|
| `HOTSPOT_NMS_MINSEP_A` | Minimum Å separation between detected hotspots when applying non-maximum suppression. Prevents duplicate centers in dense regions. |
| `R_MIN_CAVITY_A` | Minimum inscribed sphere radius (Å) required for a cavity to be considered viable. Filters out shallow surface pockets. |
| `SURFACE_SHELL__MIN_A` / `SURFACE_SHELL__MAX_A` | Inner/outer Å bounds for the surface shell used to classify near-surface voxels. |
| `SURFACE_NMS_MINSEP_A` | Non-maximum suppression distance (in voxels) when evaluating surface cavities. Larger values merge nearby openings. |
| `MAX_CENTER_DIST_A` | Å-distance threshold from the protein surface for accepting automatically detected centers. |
| `CONTACT_SHELL_A` | Thickness of the contact shell (in voxels ≈ Å) counted when evaluating pocket accessibility. |
| `HOTSPOT_BOX_ANGLE` | Minimum side length (Å) for the automatically generated search box, ensuring consistent grid volumes even for narrow cavities. |
| `MIN_SURFACE_FRAC` | Minimum fraction of grid voxels that must belong to the surface shell for a box to qualify as a surface pocket. |
| `AUTOSITES` | Target number of hotspots (grid boxes) to generate per receptor when running in automatic centers mode. |

Tweak these parameters only when you need to bias the hotspot finder—for
example, tightening `MIN_SURFACE_FRAC` to focus on buried cavities or lowering
`AUTOSITES` to restrict the number of generated docking boxes.

> **Note:** Values expressed in voxels (e.g., `SURFACE_NMS_MINSEP_A` and
> `CONTACT_SHELL_A`) can be converted to Ångström by multiplying by
> `GRID_SPACING`.

---

## Pipeline Segments & What Makes Ultidock Different

Ultidock is organized into four primary segments. Each segment has been
engineered for reliability and reproducibility compared to ad-hoc docking
scripts.

1. **Setup (`setup.py`)**
   - Idempotently creates the full directory tree (LIGANDS, MACRO_MOL, DOCKING,
     RESULTS, etc.).
   - Detects GPU availability and compiles AutoDock-GPU/AutoGrid with the
     correct compute capabilities.
   - Respects explicit CLI paths so scripted runs can reuse shared toolchains.

2. **Ligand Preparation**
   - `ligands.wget` entries are executed with robust retry logic and optional
     HTTPS upgrades (HSTS aware) unless `--skip-wget` is specified, in which case
     pre-seeded ligand archives are used as-is.
   - `extract.py` orchestrates AutoDock Vina's `vina_split` to extract, split,
     and stage ligands with deterministic filenames so downstream consumers can
     glob without guessing naming schemes.

3. **Docking (`dock_v02.py`)**
   - Per-receptor grid caching eliminates redundant AutoGrid runs even when the
     pipeline is restarted.
   - Semaphore-guarded worker pool maintains one AutoDock-GPU process per GPU
     while CPU preparation remains concurrent.
   - Metadata (grid centers, hotspots, cavity statistics) is persisted for MD
     seeding and reproducibility.

4. **I/O Hardening (`molguard`)**
   - Fixed-width float formatter (`fixedfmt.py`) ensures every number written to
     AutoGrid/AutoDock files respects the Fortran-style column widths those tools
     parse — no exponent notation, no missing decimals, no locale drift.
   - PDBQT linter and normalizer catches column-format bugs before they reach
     AutoGrid, with fail-slow error collection and a regression test suite.
   - Receptor canonicalizer produces deterministic, byte-identical files across
     machines given the same input, simplifying reproducibility audits.
   - Grid map checker validates AutoGrid output immediately after each run:
     all-zero maps, NaN/Inf energies, and missing files are caught with
     actionable error messages pointing to the `.glg` log.

5. **Analysis (`analyse_docking_results.py`)**
   - Optional stage that aggregates top poses, binding energies, and summary
     statistics. If `pandas` is unavailable the pipeline logs a warning and
     continues so production runs are never blocked by optional tooling.
   - Results are parsed directly from the automatically maintained SQLite
     database so reruns can resume and analytics scripts can attach without
     bespoke exports.

- **Single-command automation:** `run.py` orchestrates everything from toolchain
  compilation to final scoring, eliminating manual multi-step checklists.
- **Directory-first design:** explicit, user-configurable directories keep
  receptors, ligands, grids, and results isolated and reproducible.
- **Deterministic I/O:** the `molguard` layer guarantees that the same input
  always produces byte-identical PDBQT and grid files regardless of machine or
  locale, enabling reliable comparative studies.
- **Example-driven:** the `examples/` directory demonstrates full CPU and GPU
  runs, including workspace reset, staging, and pipeline invocation.
- **Resilient defaults:** built-in fallbacks for missing optional dependencies
  (e.g., pandas, wget SSL issues) keep long batches running with informative
  warnings.
- **Database-native:** every docking job streams its status into the SQLite
  results store, enabling instant post-processing without manual log parsing.
- **Future-ready:** the branch maintains alignment with planned GROMACS
  integration by preserving metadata required for MD restarts and analysis.

---

## Spotlight: Grid Boxing & Cavity Finder Algorithm

Ultidock proudly features a **high-precision grid boxing and cavity finder
algorithm** that automatically identifies docking hotspots, sizes grids to the
appropriate search volume, and surfaces cavity statistics for every receptor.
This algorithm—implemented in [`docking/dock_v02.py`](docking/dock_v02.py)—is a
cornerstone of what sets Ultidock apart. Special acknowledgement goes to the
original contributors who engineered the routine: their work enables the
pipeline to deliver reproducible, multi-site docking without manual box tuning.

- Multi-scale hotspot detection guards against missed binding pockets even on
  flexible receptors.
- Adaptive bounding boxes trim unneeded search space, speeding up AutoGrid and
  AutoDock-GPU runs while preserving accuracy.
- Persisted cavity metadata (center coordinates, occupancy metrics, and grid
  spacing) feeds directly into MD seeding and downstream analysis.

Whenever you run Ultidock, this algorithm silently prepares precise search
volumes so that the subsequent docking stages focus on the most promising
regions.

---

## Working with the Example Pipelines

Two curated examples (`gabaa-benzos` and `sert-escitalopram`) showcase the full
workflow. Each example runner performs the same steps a user would follow:

```bash
python3 examples/sert-escitalopram/example-run.py
```

What the helper (`examples/common.py`) does:

1. Calls `python3 docking/clean.py -y --all` to ensure a fresh workspace.
2. Recreates the canonical directories under `docking/`.
3. Copies the example receptor and ligands into the main pipeline directories.
4. Executes `python3 docking/run.py` with explicit path overrides.

Use these scripts as blueprints for your own automation or CI workflows.

---

## Troubleshooting

- **SSL errors while downloading ligands**
  - Corporate firewalls or strict TLS inspection can block `files.docking.org`.
    Download the required archives manually and place them in
    `docking/LIGANDS_DIR/` before running the pipeline.

- **AutoDock-GPU compilation failures**
  - Ensure CUDA 12.8+ is installed and `nvcc --version` reports the expected
    toolkit. Re-run `python3 docking/clean.py -y --all` followed by
    `python3 docking/run.py --mode gpu`.

- **`ultidock doctor` shows `[WARN] not compiled` for AutoGrid or AutoDock-GPU**
  - The source tree is present but the binaries have not been built yet.
    Run `cd docking && python setup.py` to compile them. After a successful
    build, `doctor` will report `[OK]` with the resolved binary path.

- **`ultidock pdbqt check` reports `NO_DECIMAL` or `EXPONENT` errors**
  - These indicate the PDBQT file was written by a tool that does not respect
    AutoDock column widths (e.g., some OpenBabel versions or AMBER converters).
    Run `ultidock pdbqt normalize <file> -o <fixed.pdbqt>` to reformat the
    numeric columns before docking.

- **Optional analysis skipped**
  - If you see `ModuleNotFoundError: pandas`, install it with
    `pip install pandas` and re-run the analysis stage:
    `python3 docking/analyse_docking_results.py`.

- **Out-of-disk-space errors**
  - Ligand archives can be large. Clean up with `python3 docking/clean.py -y` or
    remove unused files from `docking/LIGANDS_DIR/` and `docking/DOCKING_DIR/`.

---

## Citation & License

If you use Ultidock in academic or industrial research, please cite:

> Turgut, T. (2025). *Ultidock: A Lightweight Parallelized Docking Pipeline for
> Ligand Screening*. GitHub Repository. https://github.com/taka78/ultidock

Ultidock is released under the [MIT License](LICENSE). When applicable, please
also cite:

- Trott, O., & Olson, A. J. (2010). *AutoDock Vina: Improving the speed and
  accuracy of docking with a new scoring function, efficient optimization, and
  multithreading*. Journal of Computational Chemistry, 31(2), 455–461.
- Santos-Martins, D., et al. (2021). *Accelerating AutoDock4 with GPUs and
  Gradient-Based Local Search*. Journal of Chemical Theory and Computation,
  17(2), 1060-1073.

---

*If you find Ultidock useful, please star the repository and consider sharing
your improvements via pull requests.*
