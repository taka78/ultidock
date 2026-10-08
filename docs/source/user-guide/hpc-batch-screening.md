---
orphan: true
---

# Independent concurrent runs

One Ultidock process already schedules all ligands in its library across a
bounded CPU or GPU worker pool. Start there with
[Batch screening with Ultidock](../tutorials/hpc-screening.md). This page is
for the different case where several **independent Ultidock processes** run
at the same time. Complete a representative screen first to estimate time,
memory and disk from the [screening walkthrough](../tutorials/virtual-screening.md).

## Run folders and shared state

High-level commands (`ultidock cavity`, `blind`, `known-site`, `fpocket` and
`p2rank`) automatically create a timestamped run folder by default. Each folder contains
its own `docking/`, `analysis/` and `results/` directories; the SQLite database
is created at `results/ultidock_results.db`. You do not need to partition those
outputs by hand. Use `--output-dir` only when you want a specific run folder;
choose a different path for each concurrent job.
Bundled example scripts also create a separate workspace for each run. The
lower-level `ultidock run` command uses configured directories directly, so
provide distinct paths if you use it for concurrent jobs.

The pipeline still generates `docking/config.py` in the active application
workspace. Concurrent jobs must not share that writable configuration: give
each job its own `ULTIDOCK_HOME`. Grid generation also writes alongside the
receptor, so stage the receptor in a separate writable directory per process. Use
absolute paths and stage only the ligands intended for each process. The installed CLI
materializes workflow resources in each application's home.

## Budget concurrency across processes

Within each process, `ULTIDOCK_WORKERS` sets the maximum concurrent ligand
workers; it does not limit the total number of ligands processed. CPU runs also
use `VINA_CPU` threads per Vina process. CUDA runs distribute work across
detected NVIDIA GPU IDs using `GPU_SLOTS_PER_DEV` (default 2); OpenCL uses the
runtime's visible device.
The runner derives an OpenMP budget from host CPU count and worker count; do not
assume it automatically interprets every scheduler allocation. Bind tasks to the
allocated cores and check logs for actual worker/thread settings.

The high-level commands' separate run folders keep SQLite records from
unrelated jobs apart. Grid maps can be large; account for storage, scratch
lifetime and file-transfer cost. Copy each complete run folder and its
preparation metadata back before scratch is removed. Avoid running cleanup
against a directory still in use.

See [Batch screening with Ultidock](../tutorials/hpc-screening.md) for the
normal one-process workflow and the
[configuration reference](../reference/configuration.md) for which settings
are CLI flags versus generated Python variables.
