# Step-by-step high-throughput screening

Ultidock can screen a ligand library with one pipeline command. It prepares
receptors, selects the available docking backend, proposes binding sites,
builds grids, processes ligand archives, docks and analyzes results. Install
the [required packages and GPU runtime](../getting-started/requirements.md)
first. For a small run using supplied molecules, see
[Run the bundled examples](examples.md).

**A direct library source:** choose a chemically appropriate 3D tranche set
in [ZINC20](https://zinc20.docking.org/tranches/home/) or
[ZINC-22](https://cartblanche22.docking.org/), export **AutoDock PDBQT.gz**
with **WGET** commands, and save that command list as
`docking/ligands.wget`. Ultidock downloads and splits the selected archives
when you run `ultidock run`. The [ZINC instructions](../user-guide/start-docking.md)
show the exact file handoff and format checks. ZINC's
[3D archive](https://cache.docking.org/3D/) describes the available tranche
properties and PDBQT files.

## 1. Place the inputs

Use the active workspace printed by `ultidock doctor`. Put receptor `.pdb`,
`.mol2` or `.pdbqt` files in its `docking/MACRO_MOL_DIR/`. For ligands, use
one of these paths:

| Input | Place it here | What the pipeline does |
| --- | --- | --- |
| Prepared ligand PDBQT files or `.pdbqt.gz` archives | `docking/LIGANDS_DIR/` | Uses those files; pass `--skip-wget` for a local-only screen |
| Selected ZINC AutoDock PDBQT download list | Save the exported `wget` commands as `docking/ligands.wget` | Runs by default, downloads archives into `LIGANDS_DIR`, then extracts and splits them |

`ligands.wget` is useful for repeating a screen from the same selected
archive URLs. The bundled file contains one example download; replace it with
your selection. The current pipeline runs it
even when local ligands are present, so use `--skip-wget` to screen only local
files. Pass `--wget FILE` to select a different manifest. [Dock your own molecules](../user-guide/start-docking.md)
explains the manifest format and input directories. The pipeline checks and
prepares staged receptors; it expects ligand PDBQT chemistry to be prepared
already.

## 2. Start the screen

From an activated Ultidock environment:

```bash
ultidock run
```

If you supplied local ligands instead of a download list, run
`ultidock run --skip-wget` so the bundled example ligand is not added.

`run` executes setup, receptor preparation, the selected ligand download,
site discovery, grid generation, docking and analysis. Its default hardware
mode is `auto`: it selects CUDA when a compatible NVIDIA GPU is visible,
otherwise a visible OpenCL GPU, otherwise CPU Vina. When several NVIDIA GPUs
are detected, AutoDock-GPU jobs are distributed across them. The default is
two concurrent slots per detected device; the terminal log prints the selected
backend, GPU IDs, workers and slots. A GPU build also needs its driver and
compute toolkit; see [System requirements](../getting-started/requirements.md).

If a screen must use GPU or fail instead of falling back to CPU, run
`ultidock run --mode gpu`. If you need a reproducible comparison across
machines, record the selected backend and settings from the log. CPU Vina and
AutoDock-GPU use different search and scoring paths; compare ranked results
within a consistent backend.

The default site search is CaV-EMPS (`--grid-mode centers`). To search one
broad box instead, use `ultidock run --grid-mode blind`. Optional fpocket and
P2Rank methods have their own commands and tools; see
[Binding-site discovery](../user-guide/binding-site-discovery/index.md).

## 3. Find and review the results

The command prints the active folders and discovered ligand count. In the
workspace, `docking/MACRO_MOL_DIR/` contains prepared receptors and grids,
`docking/DOCKING_DIR/` contains 3D output poses, and
`docking/RESULTS_DIR/` contains the SQLite score database and CSV tables.
[Results and reports](../user-guide/results-reports.md) explains how to find
the scored pose corresponding to a ligand. A filtered CSV can be empty even
when docking completed; inspect the database and logs before interpreting it.

## Plan capacity for a large library

No fixed ligands-per-hour figure applies across receptors, ligands, site
counts and hardware. Use completed ligand–site jobs per wall hour from a
representative run on the machine that will screen the full library. A first
estimate is `planned ligand–site pairs / measured completed pairs per hour`,
plus grid generation and analysis time. For example, 40 completed pairs in
one hour suggests about 100 hours for 4,000 pairs on the same machine and
settings. This is a planning estimate, not a completion guarantee.

| Change | What to expect | What to check |
| --- | --- | --- |
| More NVIDIA GPUs | More AutoDock-GPU jobs can run concurrently; the runner assigns jobs across detected devices | The log's GPU IDs, utilization and available GPU memory |
| `GPU_SLOTS_PER_DEV` (default 2) | More simultaneous jobs per detected GPU | Device memory and completed pairs per hour; more slots are not always faster |
| `ULTIDOCK_WORKERS` | Sets the maximum concurrent ligand workers, not the total library size | Host CPU and memory use; on a cluster, stay within the allocated resources |
| More binding sites or larger boxes | More ligand–site work and larger grids | Site coverage, disk use and run time |
| CPU Vina threads (`--vina-cpu`) | More threads per CPU job | Worker count times Vina threads versus available cores |

For the built-in worker scheduling and an optional cluster allocation, see
[Batch screening with Ultidock](hpc-screening.md). The high-level `cavity`
command creates a separate result folder automatically. Keep the
input manifest or local ligand list, selected backend, receptor preparation
notes, site settings and the resulting poses/database with any ranked screen.
