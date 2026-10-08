---
orphan: true
---

# Advanced setup and pipeline controls

For normal use, [Dock your own molecules](../user-guide/start-docking.md) runs setup,
preparation, grids, docking and analysis together. This page explains the
separate setup command, saved configuration and stage-reuse controls for
advanced workflows. For a guided three-ligand run, use
[First Docking Run](first-docking-run.md). For a large library, use the
[screening walkthrough](../tutorials/virtual-screening.md).

## 1. Prepare the active workspace

Regular installs materialize scripts, examples, benchmarks and native build
sources under `$XDG_DATA_HOME/ultidock` (normally
`~/.local/share/ultidock`). `ULTIDOCK_HOME` can select a different writable
home; `ultidock doctor` prints the active one. Editable installs and direct
checkout commands use the repository by default. Package upgrades refresh
managed sources while preserving configuration and results. Wheels contain
native sources, not machine-specific docking executables.

From a regular install, commands work in any directory. If you use
`/usr/bin/python3 -m cli.ultidock` without installing the package, run from
the repository root and substitute that form for `ultidock` below.

## 2. Run setup

The defaults detect the available GPU or CPU backend. Run setup separately
only when you need its generated configuration before docking:

```bash
ultidock setup
ultidock doctor
```

The first setup can take time: it finds or builds AutoGrid and, in GPU mode,
AutoDock-GPU. It creates working directories, writes `docking/config.py`,
checks receptors, and runs the default `ligands.wget` list even if local PDBQT
ligands or archives exist. Add `--skip-wget` for local-only setup, or
`--wget FILE` to use a different download list. `--mode gpu` requires a detected GPU instead of
allowing CPU fallback. The lower-level checkout command is
`python3 docking/setup.py`.

An available `autogrid4` on `PATH` is reused; otherwise setup builds it from
the included sources. Regular installs link available `vina` and
`vina_split` into the managed workspace. Source checkouts may use their
bundled Linux Vina binaries. Native tools and GPU drivers remain separate from
the Python package. See [System requirements](requirements.md).

## 3. Stage only the intended inputs

The default `docking/MACRO_MOL_DIR/` holds receptor `.pdb`, `.mol2` or
`.pdbqt` files. Setup prepares each discovered receptor by filename stem;
when raw and prepared forms share a stem, it prefers the PDBQT unless forced.
The default `docking/LIGANDS_DIR/` holds ligands. Every top-level
`*.pdbqt` there is selected for docking.

Use either prepared local PDBQT ligands or put `wget`
commands for `.pdbqt.gz` archives in `docking/ligands.wget` (or pass
`--wget FILE`). The bundled manifest has an example download; see
[Dock your own molecules](../user-guide/start-docking.md) for its format and
the difference between local and downloaded ligands. The extraction stage
expands archives and uses Vina's
`vina_split`; it keeps the original `.gz` by default and removes the
temporary unsplit PDBQT so it is not docked as an extra ligand. Benchmark mode
uses compact archive behavior; `--no-keep-artifacts` selects that behavior for
an ordinary run. Validate the resulting ligand count and filenames.
Use absolute `--macro-mol-dir` and `--ligands-dir` paths for scripted runs.

See [Preparing a receptor](../user-guide/preparing-receptor.md) and
[Preparing ligands](../user-guide/preparing-ligands.md).

## 4. Run and inspect

After setup, use the saved backend and inputs:

```bash
ultidock run --skip-setup
```

Or do setup and docking together:

```bash
ultidock run
```

If the active ligand directory already contains the whole library, use
`ultidock run --skip-wget` to exclude the bundled example download.

The pipeline extracts archives, builds site grids, docks and analyzes results.
After a successful extraction, `ultidock run --skip-setup --skip-extract`
reuses the existing configuration and prepared ligands. Setup flags do not
modify the saved configuration when `--skip-setup` is present. For a controlled
change to `SITE_POLICY` or another config-only value, run setup, edit the
active `docking/config.py`, then run with `--skip-setup`.

Setup prints tool paths; docking prints ligand discovery, grid preparation,
worker launches and database insertions. In a raw pipeline run, output poses
are under `DOCKING_DIR`, receptor grids under `MACRO_MOL_DIR`, and the live
SQLite database and filtered exports under `RESULTS_DIR`.
See [Results and reports](../user-guide/results-reports.md).

For downstream visualization or scoring audits, use the notebooks in
`data-analyses/` after checking the saved poses and score associations.

## Run folders and reuse

`ultidock cavity`, `blind`, `known-site`, `fpocket` and `p2rank` create
separate run folders under `docking/RESULTS_DIR/` by default.
`--output-dir PATH` chooses one explicitly. Each keeps its run description,
docking outputs, analysis, database, CSV files and report. Pocket method runs
also keep the staged receptor, predictions, raw predictor output and grids
there; other modes keep grids in their configured `MACRO_MOL_DIR`. Existing
folders under the older top-level `runs/` location are left in place.

A second batch should have dedicated input and output directories.
`ultidock clean -y` removes build products and Python caches only;
`--results` and `--maps` explicitly remove those artifacts, while
`--all` requests a full cleanup including ligands and results. Inspect
`ultidock clean --all` before confirming deletion;
`ultidock clean -y --all` executes that full reset. Concurrent runs must use
separate `ULTIDOCK_HOME` values because setup writes shared configuration.

## Validate with MolGuard

```bash
molguard pdbqt check path/to/receptor.pdbqt
molguard receptor canonicalize receptor.pdbqt -o receptor_canon.pdbqt
molguard grids check path/to/receptor.maps.fld
molguard doctor
ultidock doctor
```

With the system-Python checkout path, use
`/usr/bin/python3 -m cli.molguard` in place of `molguard` from the
repository root. The grid checker catches all-zero maps, missing files,
nonfinite energies and atom-type mismatches.
The checkout's `python3 docking/extract.py --help` describes the separate
archive-splitting utility for filtered ligand subsets.
