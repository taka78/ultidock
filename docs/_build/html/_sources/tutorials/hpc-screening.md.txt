# HPC screening example

This Slurm skeleton runs two CPU ligand shards. Adapt the partition, memory, time,
paths and tool modules to your cluster. Before submitting an array, validate one
job on a compute node. Install Ultidock in a shared virtual environment and make
Vina, Open Babel and required native build dependencies available there.

Prepare `/project/screen/shards/0` and `/project/screen/shards/1`, each containing
only its ligand PDBQT files, plus `/project/screen/receptor.pdbqt`.

```bash
#!/bin/bash
#SBATCH --job-name=ultidock
#SBATCH --array=0-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
set -euo pipefail
source /project/venvs/ultidock/bin/activate

job_root="/project/screen/jobs/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$job_root/receptors" "$job_root/ligands"
cp /project/screen/receptor.pdbqt "$job_root/receptors/"
cp /project/screen/shards/"$SLURM_ARRAY_TASK_ID"/*.pdbqt "$job_root/ligands/"

export ULTIDOCK_HOME="$job_root/application"
export ULTIDOCK_WORKERS=4
srun --cpu-bind=cores ultidock cavity --autosites 6 --mode cpu --vina-cpu 2 \
  --skip-wget --macro-mol-dir "$job_root/receptors" \
  --ligands-dir "$job_root/ligands" --output-dir "$job_root/run"
```

Each job uses its own generated configuration, receptor grids and database. Four
workers with two Vina threads each budget eight Vina CPU threads. The runner also
sets OpenMP-related environment variables using the host CPU count; inspect logs
and enforce scheduler CPU binding rather than assuming allocation-aware detection.
The memory and time requests above are starting examples, not capacity guarantees.

The independent application homes may each need a native build on first use.
For larger arrays, pre-provision a tested tool installation and pass explicit tool
directories; do not make many jobs build into the same writable source directory.
GPU allocations require your cluster's device directives and compatible runtime;
choose `--mode gpu` only after validating device visibility inside a job.

Keep each shard's raw outputs and preparation metadata. Aggregate results only
after checking failed or missing jobs and duplicate ligand identifiers. See
[HPC / batch screening](../user-guide/hpc-batch-screening.md) for workspace rules.
