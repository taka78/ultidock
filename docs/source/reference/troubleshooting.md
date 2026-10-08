# Troubleshooting

Start with `ultidock doctor`, the printed workspace path and the first
failing subprocess log. Keep the command, actual `docking/config.py`,
inputs and tool versions before retrying. A Python environment cannot
install a system GPU driver or vendor runtime.

## Installation and native builds

| Symptom | Check and action |
| --- | --- |
| `ultidock: command not found` | Activate the venv in this shell, then `python -m pip install -e .` if needed; the repository root has no `make install` target. For the system-Python checkout route use `/usr/bin/python3 -m cli.ultidock` from the root. |
| Doctor says `[WARN] not compiled` | Native source exists but the binary is absent. Run `ultidock setup` with inputs in place to build tools for the detected backend; successful builds become `[OK]` with a path. |
| AutoGrid build or map failure | Check Autotools, `m4`, Perl and `csh`; read the first error and `.glg` log. Inspect receptor atom types, grid parameters, executable and output permissions. |
| `CL/opencl.h` or `CL/cl.h` missing | Install `ocl-icd-opencl-dev` on Debian/Ubuntu. A device runtime or working `clinfo` does not supply development headers and `libOpenCL.so`. For custom SDKs check `GPU_INCLUDE_PATH` and `GPU_LIBRARY_PATH`. |
| CUDA compilation fails | Check `nvidia-smi -L`, `nvcc --version`, driver/toolkit compatibility and the `gcc-12`/`g++-12` build tools. |
| GPU not detected or `clinfo` shows zero platforms | Check the runtime and device from the same container or compute node that runs Ultidock; an ICD loader or desktop display alone is insufficient. Use `--mode cpu` only when CPU Vina is intended. |
| fpocket fails at `src/fparams.c` with `strcpy` pointer errors | Run `bash scripts/install_pocket_tools.sh fpocket` from an updated checkout; the installer patches the pinned fpocket 4.2.3 source before building. |
| P2Rank says `Unsupported class file major version 69` | P2Rank 2.5 needs Java 17–23. Install `openjdk-21-jre-headless`; unset an incompatible `JAVA_HOME` or set `JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64`. Verify the selected Java executable. |

See [System requirements](../getting-started/requirements.md) for exact
native package and GPU checks.

## Inputs, grids and long runs

| Symptom | Check and action |
| --- | --- |
| `molguard pdbqt check` reports `NO_DECIMAL` or `EXPONENT` | A converter wrote numeric fields outside AutoDock's fixed-width form. For a ligand, run `molguard pdbqt normalize FILE -o fixed.pdbqt`; this fixes formatting, not chemistry. |
| Unknown atom type, missing coordinate/charge or nonfinite value | Return to the source molecule and preparation settings; normalization cannot create missing chemical information. |
| No HD/HS donor hydrogens | Keep automatic receptor preparation enabled and install Open Babel; inspect `.prep.json` and warnings. A partially hydrogenated receptor may require full source preparation. |
| Recovery changes heavy atoms | Do not bypass the preservation check. Compare the `.molguard-backups/` original with the source and prepare a suitable receptor again. |
| Too many ligands discovered | Every top-level `*.pdbqt` in `LIGANDS_DIR` is selected. Remove unrelated files or use a dedicated `--ligands-dir` for each batch; bundled examples create isolated workspaces. |
| AutoGrid prints S1, S2, … and pauses | It builds one map set per predicted site, sequentially. Inspect the current site's `grid.glg`; docking starts after `[autogrid] prepared N site grids`. |
| Download fails with SSL/TLS errors | A firewall or TLS inspection may block `files.docking.org`. Manually download the intended archives into `LIGANDS_DIR`, verify them, then run with `--skip-wget` to avoid the bundled download. |
| Disk runs out | Check archives, map files, raw poses and scratch with `du -sh`/`df -h`. Use a new larger workspace or `--work-root` for benchmark scratch. Archive results before any explicit clean; `ultidock clean -y` alone only clears build/cache files. |
| Another run changes settings | Give concurrent jobs separate `ULTIDOCK_HOME` and writable receptor/grid directories. High-level commands already create separate run folders and result databases. |

## Results and reports

| Symptom | Check and action |
| --- | --- |
| Analysis is skipped with `ModuleNotFoundError: pandas` | Run `python -m pip install pandas` in the **same** Python environment and rerun `python3 docking/analyse_docking_results.py` with the generated configuration. |
| CSV has no rows | Check docking failures, then inspect SQLite or make an unfiltered export as shown in [Results and reports](../user-guide/results-reports.md). Headers can remain when all poses are filtered out. |
| `docking_file` points to an input ligand | Determine whether the row comes from an older database; new best-pose metadata cannot be reconstructed for every legacy record. |
| Generic HTML report has empty tables | The report generator recognizes specific site and score files, not every SQLite database or timestamped CSV. Inspect the raw database and engine CSV. |
| Scores look good but poses do not | Open the best output PDBQT with the matching prepared receptor; check site, model, clashes, preparation and scoring assumptions. A score alone does not validate a hit. |

## Useful bug reports

Include Ultidock and engine versions, operating system, backend, command,
active workspace, first relevant error, whether setup was skipped, whether a
database was reused, and a minimal shareable input. Mention whether a fresh
bundled example reproduces the failure. Do not publish credentials or
restricted molecular datasets.
