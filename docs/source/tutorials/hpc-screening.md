# HPC screening example

This Slurm skeleton runs two ligand shards, with up to two allocated NVIDIA GPUs
per job. Adapt the GPU request syntax, partition, memory, time, paths and tool
modules to your cluster. Before submitting an array, validate one job on a
compute node. Install Ultidock in a shared virtual environment and make
AutoGrid, AutoDock-GPU build tools, Vina, Open Babel and the GPU runtime available.

Prepare `/project/screen/shards/0` and `/project/screen/shards/1`, each containing
only its ligand PDBQT files, plus `/project/screen/receptor.pdbqt`.

```bash
#!/bin/bash
#SBATCH --job-name=ultidock
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
set -euo pipefail
source /project/venvs/ultidock/bin/activate

job_root="/project/screen/jobs/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$job_root/receptors" "$job_root/ligands"
cp /project/screen/receptor.pdbqt "$job_root/receptors/"
cp /project/screen/shards/"$SLURM_ARRAY_TASK_ID"/*.pdbqt "$job_root/ligands/"

export ULTIDOCK_HOME="$job_root/application"
srun --cpu-bind=cores ultidock cavity --autosites 6 \
  --macro-mol-dir "$job_root/receptors" \
  --ligands-dir "$job_root/ligands" --output-dir "$job_root/run" --skip-wget
```

`ULTIDOCK_HOME` keeps generated configuration separate, and the staged receptor
directory keeps grid files separate. `--skip-wget` prevents the bundled
example download from entering each ligand shard. The `cavity` command creates the run's
`docking/`, `analysis/` and `results/` directories and SQLite database inside
`$job_root/run` automatically. Here, `--output-dir` chooses a predictable
location for copying results off scratch; it is optional. With two visible
NVIDIA GPUs, the default worker count is four (two slots per device), and jobs
are assigned across the GPUs. The terminal prints the detected GPU IDs and
worker count. The runner also sets OpenMP-related environment variables using
the host CPU count; inspect logs and enforce scheduler CPU binding rather than
assuming allocation-aware detection. The memory and time requests above are
starting examples, not capacity guarantees. On a CPU-only partition, remove
the GPU request; automatic mode uses Vina and sizes workers from CPU count
and `VINA_CPU`. To require GPU execution, add `--mode gpu`.

The independent application homes may each need a native build on first use.
For larger arrays, pre-provision a tested tool installation and pass explicit tool
directories; do not make many jobs build into the same writable source directory.
GPU allocations require your cluster's device directives and compatible runtime.
Check `nvidia-smi -L` inside the job before scaling up.

Keep each shard's raw outputs and preparation metadata. Aggregate results only
after checking failed or missing jobs and duplicate ligand identifiers. See
[HPC / batch screening](../user-guide/hpc-batch-screening.md) for workspace rules.
