# Batch screening with Ultidock

One Ultidock run can screen a library of many ligands. Ultidock discovers the
ligand files, queues them in its own bounded worker pool and reports progress
as they finish. The number of **workers running at once** is separate from the
total number of ligands to screen. One run handles the library without manually
splitting it into jobs or creating a directory per ligand.

## Start a library screen

Install Ultidock and the [required native tools and GPU runtime](../getting-started/requirements.md).
Run `ultidock doctor` to find the active workspace. Place the receptor in its
`docking/MACRO_MOL_DIR/`. Then either put prepared ligand PDBQT files in
`docking/LIGANDS_DIR/`, or save a selected ZINC AutoDock PDBQT download list as
`docking/ligands.wget`. [Dock your own molecules](../user-guide/start-docking.md)
explains the formats and how the download list works.

For a `ligands.wget` library, run:

```bash
ultidock cavity
```

For local ligands only, skip the bundled example download:

```bash
ultidock cavity --skip-wget
```

The command prepares the receptor, proposes CaV-EMPS sites, builds grids,
submits every discovered ligand for docking and creates a timestamped result folder. The
input folders and result path are managed by Ultidock. Repeat the same command
on a workstation or a compute node with the needed tools available.

## How the worker scheduler uses hardware

The default `auto` mode selects an available GPU backend, or CPU Vina if none
is detected. In CPU mode, Ultidock sizes its worker pool from the host CPU
count and Vina threads per worker. In CUDA mode, it distributes ligand work
across detected NVIDIA GPUs; by default it allows two concurrent slots per
GPU. The log prints the selected backend, ligand count and worker count, plus
GPU IDs in GPU mode. One run can process more ligands than there are workers or GPU slots: the
remaining ligands wait for a worker.

To choose a different maximum number of concurrent ligand workers, set
`ULTIDOCK_WORKERS` to a positive integer. For example:

```bash
ULTIDOCK_WORKERS=8 ultidock cavity --skip-wget
```

This still queues **all** ligands found in `LIGANDS_DIR`; `8` limits how
many ligand workers run at once. Pick the value for the available CPU, memory
and GPU capacity. In CUDA mode, `GPU_SLOTS_PER_DEV` also bounds simultaneous
docking work on each GPU. Increasing either setting does not guarantee a faster
screen. See [performance planning](virtual-screening.md) for a pilot-based
estimate and the [configuration reference](../reference/configuration.md) for
the controls.

## Running on a cluster

If a cluster requires a batch scheduler such as Slurm, request the CPUs, memory
and GPUs you need and run the **same Ultidock command** inside one allocation.
The cluster scheduler starts the process; Ultidock schedules the ligands within
it and distributes CUDA work across the NVIDIA GPUs visible to that process.
Check GPU visibility and native tools on the compute node. For separate
Ultidock processes running at the same time, see
[Independent concurrent runs](../user-guide/hpc-batch-screening.md).
