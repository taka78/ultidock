# HPC / batch screening

First complete one representative receptor/ligand run on the actual compute node.
Confirm tool availability, GPU visibility if requested, grid generation and result
paths before submitting an array. A login-node check cannot validate a compute
node's runtime or filesystem permissions.

## Isolate jobs

The legacy pipeline generates `docking/config.py`. Concurrent runs must not write
the same application workspace. Give each job a separate `ULTIDOCK_HOME` and
separate receptor/grid, docking and results directories. The installed CLI
materializes its workflow resources in that home. Use absolute paths and stage
only the intended ligand shard.

## Budget concurrency

`ULTIDOCK_WORKERS` sets the worker pool. CPU jobs also use `VINA_CPU` threads per
Vina process. GPU jobs use `GPU_SLOTS_PER_DEV` (default 2) and detected device IDs.
The runner derives an OpenMP budget from host CPU count and worker count; do not
assume it automatically interprets every scheduler allocation. Bind tasks to the
allocated cores and check logs for actual worker/thread settings.

Separate SQLite databases per job avoid unrelated screens interleaving. Grid maps
can be large; account for storage, scratch lifetime and file-transfer cost.
Copy complete result directories and preparation metadata back before scratch is
removed. Avoid running cleanup against a directory still in use.

See the [HPC screening tutorial](../tutorials/hpc-screening.md) for a Slurm array
skeleton and the [configuration reference](../reference/configuration.md) for
which settings are CLI flags versus generated Python variables.
