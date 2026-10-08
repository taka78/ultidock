# Configuration and runtime controls

Setup writes fully resolved Python settings to `docking/config.py` in the
active application workspace. This file controls subsequent legacy pipeline
modules; it is not a YAML input schema. High-level mode wrappers also save a
`run_config.yaml` description in each run folder, but that description does
not replace the generated Python configuration. Find the workspace with
`ultidock doctor` and archive the actual `config.py` used for a publication.

## Setup and run flags

Options for `ultidock setup` and `ultidock run` are forwarded to
`docking/setup.py` and `docking/run.py`. Lowercase and legacy uppercase
directory spellings are accepted. All are optional; paths default under the
active workspace's `docking/` directory.

| Option | Default or effect |
| --- | --- |
| `--mode {auto,gpu,cuda,opencl,cpu}` | `auto` may fall back to CPU; `gpu` requires detection; other modes choose explicitly |
| `--ligands-dir` / `--LIGANDS_DIR` | Prepared ligand inputs |
| `--macro-mol-dir` / `--MACRO_MOL_DIR` | Receptor inputs, sites and grids |
| `--docking-dir`, `--analysis-dir`, `--results-dir` | Separate outputs, analysis and SQLite/export locations |
| `--vina-dir`, `--autodock-gpu-dir` | Native tool directories |
| `--wget FILE`, `--skip-wget` | Use a different download manifest or disable downloads; the bundled manifest runs by default even when local ligands exist |
| `--benchmark` | Persist benchmark-oriented settings and compact artifact behavior |
| `--grid-mode {centers,ligand,residues,blind}` | `centers` by default |
| `--grid-spacing` | 0.375 Å |
| `--grid-margin` | 5.0 Å |
| `--grid-cap` | 150.0 Å per blind-box axis |
| `--autosites` | Request 6 automatic sites |
| `--centers-tsv`, `--ref-ligand-pdb` | Site table or ligand-centered reference |
| `--vina-cpu`, `--vina-exhaustiveness`, `--vina-num-modes` | 2 threads, 8 exhaustiveness, 9 modes |
| `--vina-seed` | Unset unless supplied |
| `--receptor-prep-mode {auto,off}` | `auto` preparation/canonicalization |
| `--receptor-prepare-command TEMPLATE` | Converter with `{input}`, `{output}`, optional `{seed}` |
| `--receptor-prep-seed`, `--force-receptor-prep` | Seed 42; regenerate same-stem prepared output when requested |
| `--skip-profile` | Preserve geometry `.config.toml` sidecars |
| `--skip-setup`, `--skip-extract`, `--skip-docking`, `--skip-analysis` | Skip the named pipeline stage on `run` |
| `--keep-artifacts` / `--no-keep-artifacts` | Keep or discard original ligand archives |

`known-site` adds `--center` and `--box-size`; `cavity` adds
`--autosites`; pocket modes add `--receptor`, `--tool`,
`--grid-spacing`, and box options. Their `--output-dir` groups outputs
but does not move separately supplied input directories. Use absolute paths
for independent or scripted runs. See [CLI reference](cli.md).

## Directory and engine variables

| Variable | Meaning |
| --- | --- |
| `LIGANDS_DIR` | Staged archives and prepared ligand PDBQT files |
| `DOCKING_DIR` | Vina/AutoDock-GPU raw poses and logs |
| `ANALYSIS_DIR` | Intermediate scoring and aggregation files |
| `VINA_DIR` | `vina` and `vina_split` location |
| `AUTODOCK_GPU_DIR` | AutoDock-GPU and AutoGrid build/install directory |
| `MACRO_MOL_DIR` | Receptors, generated grids and per-site artifacts |
| `RESULTS_DIR` | SQLite database and final CSV/JSON exports |
| `DB_PATH` | Absolute `ultidock_results.db` path |
| `GPU_TYPE` | Selected `CPU`, `CUDA` or OpenCL backend value |
| `NUMWI` | AutoDock-GPU work-item build setting; larger values may better saturate a large GPU while smaller values may help memory-constrained devices, but changing it requires the corresponding native binary and measurement |
| `AUTO_GRID_BIN` | Resolved `autogrid4` path; adjust for a tested prebuilt installation |
| `VINA_CPU`, `VINA_EXHAUSTIVENESS`, `VINA_NUM_MODES`, `VINA_SEED` | Vina thread, sampling, output and seed settings |
| `RECEPTOR_PREP_MODE`, `RECEPTOR_PREP_COMMAND`, `RECEPTOR_PREP_SEED`, `RECEPTOR_PREP_FORCE` | Shared receptor preparation controls |

## Search-region variables

| Variable | Meaning |
| --- | --- |
| `GRID_MODE` | `centers`, `ligand`, `residues` or `blind` strategy |
| `SITE_POLICY` | `receptor_search`, `exhaustive_search`, `internal`, `surface` or `hybrid` |
| `GRID_SPACING` | Grid resolution in Å; finer spacing increases grid work |
| `GRID_MARGIN` | Å padding around hotspot-derived boxes |
| `GRID_CAP` | Maximum blind-box side length in Å |
| `CENTERS_TSV` | Receptor/site centers table, generated when needed |
| `REF_LIGAND_PDB` | Reference ligand for `GRID_MODE="ligand"` |
| `AUTOSITES` | Requested sites per receptor, not guaranteed actual count |

`SITE_POLICY` and some advanced dials lack public setup flags. Generate
configuration with `ultidock setup`, edit its `docking/config.py`, then use
`ultidock run --skip-setup`. A new setup would overwrite those edits;
`--skip-setup` also means newly supplied setup flags cannot change the
existing configuration.

## CaV-EMPS tuning variables

These control the receptor-derived site finder described in
[CaV-EMPS methodology](../scientific-background/cav-emps-methodology.md).
Change them only with a declared evaluation plan.

| Variable | Meaning |
| --- | --- |
| `HOTSPOT_NMS_MINSEP_A` | Minimum separation between candidate centers, clamped relative to box size |
| `R_MIN_CAVITY_A` | Cavity inscribed-radius threshold; `None` enables receptor-adaptive estimation |
| `ADAPTIVE_R_MIN_*` | EDT threshold estimation and clamping controls |
| `MAPS_POCKET_MAX_A` | Upper EDT shell bound for map-driven surface pockets |
| `HOTSPOT_NMS_BOX_FRACTION`, `HOTSPOT_NMS_MIN_A`, `HOTSPOT_NMS_MAX_A` | Box-relative inter-site separation bounds |
| `SURFACE_SHELL__MIN_A`, `SURFACE_SHELL__MAX_A` | Near-surface shell bounds |
| `SURFACE_NMS_MINSEP_A` | Separation floor among surface candidates |
| `MAX_CENTER_DIST_A` | Maximum center distance from protein surface |
| `CONTACT_SHELL_A` | Contact shell thickness for pocket accessibility |
| `HOTSPOT_BOX_ANGLE` | Minimum automatically generated box side in Å |
| `MIN_SURFACE_FRAC` | Minimum near-surface voxel fraction for surface pockets |

## Concurrency and workspace isolation

`ULTIDOCK_HOME` chooses the managed application workspace.
`ULTIDOCK_WORKERS` overrides the maximum number of concurrent ligand workers;
all discovered ligands are still queued for processing. CPU mode's default uses host
CPU count divided by `VINA_CPU`; in a scheduler allocation, set workers
explicitly so workers × Vina threads fits allocated cores. GPU mode uses
`GPU_SLOTS_PER_DEV` (default 2) and detected device IDs. The runner also
derives an OpenMP budget from host CPU count; check logs and bind allocated
cores rather than assuming it knows your scheduler allocation.

Independent concurrent Ultidock processes need separate generated configuration
and writable receptor/grid directories. High-level commands automatically create separate
run folders and databases; the lower-level `ultidock run` command uses its
configured output paths directly. See
[Independent concurrent runs](../user-guide/hpc-batch-screening.md) and
[performance planning](../tutorials/virtual-screening.md).
