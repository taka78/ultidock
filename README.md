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
5. [Deterministic Receptor Input Handling](#deterministic-receptor-input-handling)
6. [Configuration Reference](#configuration-reference)
7. [Pipeline Segments & What Makes Ultidock Different](#pipeline-segments--what-makes-ultidock-different)
8. [Spotlight: Grid Boxing & Cavity Finder Algorithm](#spotlight-grid-boxing--cavity-finder-algorithm)
9. [Benchmarking](#benchmarking)
10. [Working with the Example Pipelines](#working-with-the-example-pipelines)
11. [Troubleshooting](#troubleshooting)
12. [Citation & License](#citation--license)

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
  automake autoconf libtool m4 perl pkg-config \
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

Install all required packages plus the `ultidock` and `molguard` CLIs in one step:

```bash
pip install -r requirements.txt
pip install -e .
```

This installs numpy, scipy, psutil, matplotlib, pandas, the `ultidock` workflow
CLI, and the `molguard` deterministic I/O CLI.
`pandas` and `matplotlib` are used only by the post-run analysis stage.

If you want Ultidock to prepare raw receptor `.pdb` files automatically, install
at least one receptor conversion backend:

```bash
pip install meeko
# or install Open Babel through your system package manager
```

Existing receptor `.pdbqt` files do not need Meeko/Open Babel; they are only
canonicalized by `molguard`.

---

## Repository Layout

```
ultidock/
├─ pyproject.toml            # Package manifest; defines `ultidock` and `molguard` CLIs
├─ requirements.txt          # Runtime dependencies (numpy, scipy, click, ...)
├─ requirements-dev.txt      # Development/test dependencies
├─ SETUP.md                  # Step-by-step new-user guide
├─ cli/                      # Repo-root CLI entry points
│  ├─ ultidock.py            # Workflow, benchmark, and example commands
│  └─ molguard.py            # File/receptor/grid safety commands
├─ molguard/                 # I/O hardening + validation layer
│  ├─ io/
│  │  ├─ fixedfmt.py         # Fixed-width, locale-safe float formatter
│  │  ├─ pdbqt.py            # PDBQT linter, normalizer, receptor canonicalizer
│  │  └─ receptor_prep.py    # Shared receptor prep for pipeline + benchmarks
│  ├─ grids/
│  │  └─ check.py            # AutoGrid .fld / .map sanity checker
│  ├─ cli.py                 # Compatibility shim for older imports
│  └─ tests/                 # Unit and regression tests (pytest)
├─ docking/
│  ├─ run.py                 # Main entry point for the entire pipeline
│  ├─ setup.py               # Idempotent environment + dependency setup
│  ├─ dock_v02.py            # AutoDock-GPU / AutoGrid orchestration
│  ├─ make_grids.py          # Receptor-only site finder and grid generator
│  ├─ profile_receptors.py   # Optional per-receptor sidecar profiler
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
├─ benchmarks/               # DUD-E download, preparation, site recovery, and docking benchmarks
├─ examples/                 # Self-contained example runners
└─ data-analyses/, results/  # Optional downstream notebooks & exports
```

The workflow assumes you copy receptor `.pdb` or `.pdbqt` files into
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
   - Copy your receptor(s) to `docking/MACRO_MOL_DIR/`. Ultidock scans files by
     extension: each `*.pdb` is sanitized and converted to `<input-stem>.pdbqt`,
     and each `*.pdbqt` is only canonicalized. Generated receptor folders use
     the same discovered stem, matching the `dock_v02.py`/`make_grids.py` flow.
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

When a CLI command hits a common failure path, it prints a short pointer to the
most relevant README section instead of burying the terminal in long guidance.

| Command | Purpose |
|---------|---------|
| `python3 docking/run.py [options]` | Primary entry point. Validates the environment, runs setup, downloads ligands, launches docking, and triggers analysis. |
| `python3 docking/setup.py [options]` | Runs the setup stage only (directory creation, AutoDock-GPU/AutoGrid/Vina checks). All CLI flags mirror `run.py`. |
| `python3 docking/dock_v02.py [options]` | Executes the docking stage against prepared ligands and receptors. Used internally by `run.py`. |
| `python3 docking/extract.py` | Wrapper around AutoDock Vina's `vina_split` for splitting ligand archives and optional filtering. |
| `python3 docking/clean.py -y --all` | Removes compiled binaries, cached grids, downloads, and generated configs. Use before starting a fresh run. |
| `python3 docking/profile_receptors.py` | Generate optional per-receptor `.config.toml` sidecars from receptor geometry. |
| `python3 benchmarks/cavity_recovery_benchmark.py` | Evaluate receptor-only site finding against co-crystallized ligand centers. |
| `python3 benchmarks/download_dude.py` | Download DUD-E receptor, crystal ligand, active, and decoy files. |
| `ultidock run [options]` | Run the full docking pipeline from anywhere in the repo (no need to `cd docking/`). Forwards all flags to `docking/run.py`. |
| `ultidock setup [options]` | Run the setup stage through the workflow CLI. |
| `ultidock clean [-y] [--all]` | Reset compiled binaries and outputs. Forwards all flags to `docking/clean.py`. |
| `ultidock benchmark cavity-recovery [options]` | CLI wrapper for the receptor-only site-recovery benchmark. |
| `ultidock benchmark download-dude [options]` | Download DUD-E receptor, crystal ligand, active, and decoy files. |
| `ultidock example list` / `ultidock example run <name>` | Discover and run bundled example pipelines. |
| `ultidock doctor` | Print tool locations and versions. Distinguishes between binaries not compiled yet (source present) and not found at all. |
| `molguard pdbqt check <file>` | Lint a receptor or ligand PDBQT for AutoDock column-format issues (exponent notation, missing decimals, bad atom types). |
| `molguard pdbqt normalize <file> -o <out>` | Rewrite all numeric columns in a ligand PDBQT through the fixed-width formatter. Torsion tree is left untouched. |
| `molguard receptor canonicalize <file> -o <out>` | Sort, renumber, and reformat a receptor PDBQT deterministically. Returns a SHA-256 digest for reproducibility checks. |
| `molguard receptor prepare <file> -o <out>` | Run the shared receptor-prep path. `.pdbqt` is canonicalized; `.pdb` is sanitized, converted, then canonicalized. |
| `molguard grids check <maps.fld>` | Validate AutoGrid output: checks for all-zero maps, NaN/Inf energies, missing files, and atom-type mismatches. |
| `molguard doctor` | Print MolGuard version and optional receptor-conversion backend availability. |

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
| `--skip-profile` | Preserve existing per-receptor `.config.toml` sidecars and skip automatic profiling. |
| `--receptor-prep-mode {auto,off}` | Enable or disable shared receptor preparation/canonicalization. |
| `--receptor-prepare-command TEMPLATE` | Override the receptor conversion command. Use `{input}`, `{output}`, and optionally `{seed}`. |
| `--force-receptor-prep` | Regenerate prepared receptor `.pdbqt` outputs even if sibling outputs already exist. |

All flags are optional; defaults point to directories within `docking/`.

---

## Deterministic Receptor Input Handling

Ultidock does not require receptors to use a fixed filename such as
`receptor.pdb` or `gpcr_beta.pdbqt`. The pipeline scans `MACRO_MOL_DIR` by file
extension and treats the discovered filename stem as the receptor identity used
for generated receptor files, site folders, grid caches, and downstream result
names.

The shared receptor-prep path lives in `molguard.io.receptor_prep` and is used
by both regular pipeline runs and benchmark scripts:

- `*.pdbqt` inputs are not reconverted. They are canonicalized through
  `molguard` so atom ordering, numbering, fixed-width numeric fields, and the
  resulting SHA-256 digest are deterministic.
- `*.pdb` inputs are sanitized for receptor-conversion tools, converted to
  `<input-stem>.pdbqt`, then canonicalized through the same PDBQT path.
- Multiple receptor files can live in `MACRO_MOL_DIR`; each receptor keeps its
  own discovered stem. If both raw and prepared forms exist for the same stem,
  the existing `.pdbqt` is preferred unless receptor prep is forced.
- Drastic rescue actions, such as deleting residues after a Meeko excess-bond
  failure, print a loud warning because they change the receptor model and must
  be reported with benchmark or docking results.

This structure keeps the benchmark and production pipeline scientifically
aligned: receptor preparation is not a special benchmark-only script, and a
researcher can use normal receptor filenames without manually renaming files to
match Ultidock internals.

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
| `SITE_POLICY` | Automatic site-family policy: `receptor_search`, `exhaustive_search`, `internal`, `surface`, or `hybrid`. |
| `GRID_SPACING` | Ångström spacing between grid points. Smaller values yield finer resolution at the cost of longer AutoGrid runtimes. |
| `GRID_MARGIN` | Extra Ångström padding applied to each hotspot-derived grid to ensure the box fully encloses the binding site. |
| `GRID_CAP` | Maximum Å-length per axis when running in blind mode to prevent runaway grid sizes. |
| `CENTERS_TSV` | Optional path to a precomputed `centers.tsv`. Leave as `None` to let Ultidock regenerate hotspot centers automatically. |
| `REF_LIGAND_PDB` | Reference ligand file used when `GRID_MODE="ligand"` to seed the search box from a co-crystal pose. |
| `RECEPTOR_PREP_MODE` | `auto` prepares/canonicalizes receptor inputs in `MACRO_MOL_DIR`; `off` leaves files untouched. |
| `RECEPTOR_PREP_COMMAND` | Optional receptor conversion command template for `.pdb`/`.mol2` inputs. |
| `RECEPTOR_PREP_SEED` | Seed forwarded to external receptor-conversion commands. |
| `RECEPTOR_PREP_FORCE` | Regenerate converted receptor outputs when `True`. |

### Grid Boxing & Cavity Finder Dials

These parameters feed the hotspot detection routine acknowledged in the
[Spotlight](#spotlight-grid-boxing--cavity-finder-algorithm) section.

| Variable | Description |
|----------|-------------|
| `HOTSPOT_NMS_MINSEP_A` | Minimum Å separation between detected hotspots when applying non-maximum suppression. The effective value is clamped relative to box size. |
| `R_MIN_CAVITY_A` | Minimum inscribed sphere radius (Å) required for a cavity to be considered viable. `None` enables adaptive receptor-specific estimation. |
| `ADAPTIVE_R_MIN_*` | Parameters controlling receptor-specific EDT threshold estimation and clamping. |
| `MAPS_POCKET_MAX_A` | Upper EDT shell bound used by map-driven surface pocket detection. |
| `HOTSPOT_NMS_BOX_FRACTION`, `HOTSPOT_NMS_MIN_A`, `HOTSPOT_NMS_MAX_A` | Bounds used to clamp inter-site separation from the requested docking-box side length. |
| `SURFACE_SHELL__MIN_A` / `SURFACE_SHELL__MAX_A` | Inner/outer Å bounds for the surface shell used to classify near-surface voxels. |
| `SURFACE_NMS_MINSEP_A` | Å separation floor when evaluating surface cavities. Larger values merge nearby openings. |
| `MAX_CENTER_DIST_A` | Å-distance threshold from the protein surface for accepting automatically detected centers. |
| `CONTACT_SHELL_A` | Thickness of the contact shell counted when evaluating pocket accessibility. |
| `HOTSPOT_BOX_ANGLE` | Minimum side length (Å) for the automatically generated search box, ensuring consistent grid volumes even for narrow cavities. |
| `MIN_SURFACE_FRAC` | Minimum fraction of grid voxels that must belong to the surface shell for a box to qualify as a surface pocket. |
| `AUTOSITES` | Target number of hotspots (grid boxes) to generate per receptor when running in automatic centers mode. |

Tweak these parameters only when you need to bias the hotspot finder—for
example, tightening `MIN_SURFACE_FRAC` to focus on buried cavities or lowering
`AUTOSITES` to restrict the number of generated docking boxes.

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
   - Runs shared receptor preparation through `molguard`: existing `.pdbqt`
     receptors are canonicalized, while raw `.pdb` receptors are sanitized,
     converted, and then canonicalized.
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
   - Site folders are named from the discovered receptor stem; site IDs such as
     `S1` and `S2` are output identifiers only, not quality labels.
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
   - Shared receptor preparation produces deterministic receptor PDBQT files and
     prints loud warnings for drastic rescue actions, such as deleting residues
     after a Meeko excess-bond failure.
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

Ultidock's automatic site finder is one of the pipeline's central features. Its
goal is to approximate the binding-site center that a researcher might otherwise
take from a co-crystal ligand, using only receptor-derived information. In
practical screening work, that co-crystallized ligand and its centroid are often
not available, so Ultidock does not require a known crystal center before
docking.

The finder is implemented in [`docking/make_grids.py`](docking/make_grids.py)
and orchestrated by [`docking/dock_v02.py`](docking/dock_v02.py). It analyzes
receptor geometry and receptor-derived AutoGrid signals to propose a compact set
of likely docking boxes. In benchmarks, the crystal ligand center is used only
after site generation as an external recovery reference, not as an input to the
site finder.

The default `receptor_search` policy combines complementary receptor-derived
signals:

- **Internal geometry:** receptor atoms are rasterized onto the AutoGrid lattice,
  and an Euclidean distance transform identifies enclosed cavities and channels.
- **Surface/map favorability:** favorable AutoGrid interaction maps are clipped
  to physically useful negative values, smoothed over a ligand-sized region, and
  masked to receptor-proximal pocket shells. This avoids selecting a single
  favorable surface voxel as if it were a pocket center.
- **Hybrid portfolio assembly:** internal, surface, and consensus candidates are
  de-duplicated by Å-scale non-maximum suppression and written as docking boxes.

Adaptive clamping is used in two places. First, the cavity radius threshold can
be inferred from each receptor's own EDT peak distribution instead of forcing a
single global value. Second, inter-site separation is clamped relative to the
docking-box side length so the generated boxes are distinct but not needlessly
sparse.

Site IDs (`S1`, `S2`, etc.) are only output identifiers. They are not scores,
priorities, or quality labels. Downstream evaluation should use distance or
docking metrics over all generated sites rather than interpreting `S1` as the
"best" site.

---

## Benchmarking

Benchmark scripts live in `benchmarks/`. The repository tracks the scripts, but
downloaded DUD-E datasets and generated benchmark outputs are ignored by Git and
should remain local.

### Download DUD-E Targets

Download receptor, crystal ligand, active, and decoy files for selected targets:

```bash
ultidock benchmark download-dude \
  --targets ace,bace1,braf,cdk2,cxcr4,drd3,egfr,esr1,gcr,hdac2,hivpr,pde5a,pparg,src,vgfr2 \
  --dataset-root benchmarks/datasets
```

### Receptor Site-Recovery Benchmark

Evaluate whether the receptor-only site finder recovers the experimentally
observed co-crystal pocket neighborhood:

```bash
ultidock benchmark cavity-recovery \
  --dataset-root benchmarks/datasets \
  --targets ace,bace1,braf,cdk2,cxcr4,drd3,egfr,esr1,gcr,hdac2,hivpr,pde5a,pparg,src,vgfr2 \
  --autosites 6 \
  --site-policy receptor_search \
  --force \
  --output-dir benchmarks/results/cavity_recovery
```

The benchmark measures the distance from each generated site center to the
centroid of the co-crystallized ligand. The ligand centroid is an external
reference for the experimentally observed bound pose; it is not used during site
generation.

Important output files:

- `summary.csv`: target-level closest-site distances and success flags
- `sites.csv`: per-site distances for every generated docking box
- per-target `centers.tsv`: the generated docking boxes used for evaluation

### Docking/Enrichment Benchmarks

`benchmarks/dude_docking_benchmark.py` and
`benchmarks/run_full_dude_benchmark.py` run full active/decoy docking arms and
collect virtual-screening metrics such as ROC-AUC, EF1, EF5, BEDROC, and logAUC.
These runs are much more expensive than the cavity-recovery benchmark and should
be treated as final validation runs rather than quick smoke tests.

---

## Working with the Example Pipelines

Two curated examples (`gabaa-benzos` and `sert-escitalopram`) showcase the full
workflow. Each example runner performs the same steps a user would follow:

```bash
ultidock example list
ultidock example run sert-escitalopram
```

What the helper (`examples/common.py`) does:

1. Calls `python3 docking/clean.py -y --all` to ensure a fresh workspace.
2. Recreates the canonical directories under `docking/`.
3. Copies the example receptor and ligands into the main pipeline directories.
4. Executes the same pipeline path as `ultidock run`, with explicit path overrides.

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

- **`molguard pdbqt check` reports `NO_DECIMAL` or `EXPONENT` errors**
  - These indicate the PDBQT file was written by a tool that does not respect
    AutoDock column widths (e.g., some OpenBabel versions or AMBER converters).
    Run `molguard pdbqt normalize <file> -o <fixed.pdbqt>` to reformat the
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
