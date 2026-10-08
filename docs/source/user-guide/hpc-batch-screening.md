# HPC / batch screening

Ultidock can screen many receptors against a ligand library in **one run**.
Each prepared ligand is assigned to a bounded worker; that worker visits every
prepared receptor and each of its proposed binding sites. Ultidock manages the
CPU or GPU worker pool and writes scored docking records from the run into one
SQLite database. You can use the same workflow on a workstation or an HPC
compute node.

## 1. Check the compute environment

Install Ultidock and its native programs using [Installation](../getting-started/installation.md)
and [System requirements](../getting-started/requirements.md). For GPU docking,
the driver and CUDA or OpenCL runtime must be visible in the environment where
the screen runs. Run `ultidock doctor` there to see the active workspace and
tool status. A small complete [first docking run](../getting-started/first-docking-run.md)
checks more than an installation or login-node check can.

## 2. Place receptors and ligands in the workspace

Use the workspace printed by `ultidock doctor`. Put **all intended receptor
files** directly in its `docking/MACRO_MOL_DIR/`; `.pdb`, `.mol2` and `.pdbqt`
are accepted. Give each receptor a distinct filename stem, such as
`kinase_a.pdb` and `kinase_b.pdb`. Setup prepares raw structures to PDBQT,
checks the inputs, and creates per-receptor site and grid files. Every
top-level prepared receptor PDBQT in that directory is included in the run.
See [Automatic receptor preparation](preparing-receptor.md) for format and
identity details.

Choose **one** ligand source:

| Source | Put it here | Run command |
| --- | --- | --- |
| Prepared local PDBQT ligands | `docking/LIGANDS_DIR/` | `ultidock cavity --skip-wget` |
| ZINC AutoDock PDBQT archive list | Save exported `wget` commands as `docking/ligands.wget` | `ultidock cavity` |

`ligands.wget` is a text list of download commands for selected ligand
archives. Ultidock downloads and splits those archives before docking. The
bundled file contains an example download, so `--skip-wget` keeps a local-only
screen limited to your staged ligands. [Dock your own molecules](start-docking.md)
shows how to select a ZINC library and write the manifest. For local inputs,
stage only the intended ligand `.pdbqt` files directly in `LIGANDS_DIR`;
[Ligand inputs](preparing-ligands.md) covers their chemical preparation.

From an editable checkout root, this stages a local library in the standard
folders; for a regular install, use the workspace path printed by `doctor`:

```bash
mkdir -p docking/MACRO_MOL_DIR docking/LIGANDS_DIR
cp /path/to/receptors/*.pdb docking/MACRO_MOL_DIR/
cp /path/to/ligands/*.pdbqt docking/LIGANDS_DIR/
ultidock cavity --skip-wget
```

Use the extensions and paths matching your actual files. For a ZINC manifest,
put its `wget` lines in `docking/ligands.wget` and run `ultidock cavity`.

## 3. Let Ultidock schedule the screen

The default `auto` mode chooses an available GPU backend, or CPU Vina when no
GPU is detected. CUDA work is distributed across detected NVIDIA GPUs. One
run uses a single selected backend; it does not mix Vina and AutoDock-GPU
scores. The log prints the backend, discovered ligand count and worker count,
plus detected NVIDIA device IDs in CUDA mode.

`ULTIDOCK_WORKERS` controls **concurrent ligand workers**, not the number of
ligands screened. The default CPU pool is sized from host CPU count and Vina
threads; the default CUDA pool uses two slots per detected NVIDIA GPU. If you
need a different concurrent worker count, set a positive value, for example:

```bash
ULTIDOCK_WORKERS=8 ultidock cavity --skip-wget
```

This still queues every discovered ligand. Each worker processes that ligand
against all receptors and their sites. The progress counter counts completed
**ligand files**, not completed receptor–ligand–site combinations. In CUDA mode,
`GPU_SLOTS_PER_DEV` also bounds active docking work on each GPU. See
[worker scheduling and cluster allocation](../tutorials/hpc-screening.md)
and the [configuration reference](../reference/configuration.md).

## 4. Find the shared SQLite results

`ultidock cavity` creates a timestamped folder under the workspace's
`docking/RESULTS_DIR/`. Inside that folder, `results/ultidock_results.db`
contains the successfully parsed docking records for **all receptors, ligands
and sites in this run**. Raw engine poses are in its `docking/` folder;
receptor sites and grids are written under `MACRO_MOL_DIR/`. Results are saved
as they are parsed, before the final analysis export. The `ligand_name`
contains the receptor stem, `binding_site` records the site, and
`docking_file` points to the pose output. Failed docking attempts do not
become scored rows; check the log for errors. Separate `cavity` invocations
create separate run folders and SQLite databases. See [Results & reports](results-reports.md)
for exports and pose inspection.

## 5. Size a large screen

The docking workload scales with the number of ligands times the **sum of
sites across all receptors**. For example, 200 receptors and 3,000 ligands
mean 600,000 receptor–ligand combinations with one site each; six sites each
would mean up to 3.6 million docking attempts. Actual sites, successful
records, time and storage vary by input and hardware. Run a representative
[screening pilot](../tutorials/virtual-screening.md)
to measure throughput and grid/output size before a large library.

On an HPC system, request CPUs, memory and any GPUs from the cluster scheduler,
then run the same Ultidock command in the allocation. Ultidock schedules the
ligands within that process. The worker pool is based on host CPU count or
visible GPUs, so check the logged settings against the allocated resources.
Grid maps and raw poses can consume substantial scratch space; retain the
complete run folder and receptor preparation metadata with the database.

### Multiple independent Ultidock processes

You only need separate application workspaces if you deliberately start
**multiple Ultidock processes at the same time**. Setup writes the shared
`docking/config.py`, and grid generation writes beside the receptors. Give
independent processes distinct `ULTIDOCK_HOME` values and writable receptor
directories so their generated files cannot collide. Each high-level command
creates its own result folder and SQLite database. This is a separate
deployment choice from screening many receptors and ligands in one run.
