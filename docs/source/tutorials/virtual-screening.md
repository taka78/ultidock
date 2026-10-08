# Step-by-step high-throughput screening

This walkthrough takes a prepared receptor and ligand library through a small
pilot, a measured capacity estimate, a full screen and result checks. Commands
assume Linux, an activated Ultidock environment and **prepared ligand PDBQT**
files. Complete [Installation and screening preparation](../getting-started/installation.md)
first. The bundled [first docking run](../getting-started/first-docking-run.md)
is a useful three-ligand check before using your own data.

## 1. Define the screen

Choose the receptor structure, chemical-state preparation protocol, site method
and backend before comparing results. CPU mode uses Vina; GPU mode uses
AutoDock-GPU, which has different search and scoring behavior. Use one backend
and fixed settings across a single ranked screen. Keep source structures and a
table mapping each PDBQT filename to its compound, stereochemistry, protonation
state and preparation settings. Distinct states need distinct filenames.

For this walkthrough, `cavity --autosites 1` proposes one CaV-EMPS site without
requiring a manually supplied center. Inspect that site before scaling. If you
have a validated known-site box, use
`ultidock known-site --center X,Y,Z --box-size SIDE` with the same input and
mode flags instead. The
[unknown-site guide](unknown-binding-site.md) describes alternative predictors.

## 2. Stage and validate inputs

Set these paths to your own **absolute** receptor file and prepared ligand
directory. The output root is created separately, so no previous ligand files
are silently included. Keep only intended `.pdbqt` files in `LIBRARY`:

```bash
export RECEPTOR=/absolute/path/to/receptor.pdbqt
export LIBRARY=/absolute/path/to/prepared-ligands
export SCREEN_ROOT="$PWD/ultidock-screen-$(date +%Y%m%d-%H%M%S)"
mkdir -p "$SCREEN_ROOT/pilot/receptors" "$SCREEN_ROOT/pilot/ligands" \
  "$SCREEN_ROOT/full/receptors"
cp "$RECEPTOR" "$SCREEN_ROOT/pilot/receptors/"
cp "$RECEPTOR" "$SCREEN_ROOT/full/receptors/"
molguard pdbqt check "$RECEPTOR"
find "$LIBRARY" -maxdepth 1 -type f -name '*.pdbqt' | wc -l
```

The count is the number of ligand files that the full run should discover.
If your receptor is a raw PDB or MOL2, prepare it first as described in
[Preparing a receptor](../user-guide/preparing-receptor.md). Validate prepared
ligands with MolGuard before a large run:

```bash
for ligand in "$LIBRARY"/*.pdbqt; do
  molguard pdbqt check "$ligand" || exit 1
done
```

Resolve failed inputs using the original chemistry and
[ligand preparation notes](../user-guide/preparing-ligands.md). Format
normalization cannot choose protonation, stereochemistry or charges for you.

Select 10–50 representative ligands for the pilot, including flexible or unusual
molecules. The command below copies the first 20 sorted filenames as a quick
mechanical smoke test; replace or supplement them with representative molecules
before using pilot timing to size the full screen:

```bash
find "$LIBRARY" -maxdepth 1 -type f -name '*.pdbqt' -print0 |
  sort -z | head -z -n 20 |
  xargs -0 -r cp -t "$SCREEN_ROOT/pilot/ligands"
find "$SCREEN_ROOT/pilot/ligands" -maxdepth 1 -type f -name '*.pdbqt' | wc -l
```

The `sort -z`, `head -z` and `cp -t` options above are GNU/Linux tools. Verify
that both counts are nonzero and that the receptor folder contains only the
intended receptor. Do not put source SDF files or unrelated PDBQTs in either
ligand directory.

## 3. Preview the pilot

Use an explicit backend. This CPU example keeps Vina's sampling settings and
site count visible. `--skip-wget` prevents the download manifest from adding
ligands. `--output-dir` groups run outputs; it does not relocate the input
directories:

```bash
ultidock doctor
ultidock cavity --autosites 1 --mode cpu --vina-cpu 2 \
  --vina-exhaustiveness 8 --vina-seed 42 --skip-wget \
  --macro-mol-dir "$SCREEN_ROOT/pilot/receptors" \
  --ligands-dir "$SCREEN_ROOT/pilot/ligands" \
  --output-dir "$SCREEN_ROOT/pilot/run" --dry-run
```

Check the printed command and `run_config.yaml`. This dry run writes the run
directory/configuration but does not execute setup, grid generation or docking.
It does not establish that the molecules are chemically correct or that the
native engines will run.

## 4. Run and inspect the pilot

Repeat the same command without `--dry-run`. Bash `time` reports elapsed wall
time; record it together with CPU/GPU model, RAM, backend, tool versions, ligand
count, site count, box size and settings:

```bash
time ultidock cavity --autosites 1 --mode cpu --vina-cpu 2 \
  --vina-exhaustiveness 8 --vina-seed 42 --skip-wget \
  --macro-mol-dir "$SCREEN_ROOT/pilot/receptors" \
  --ligands-dir "$SCREEN_ROOT/pilot/ligands" \
  --output-dir "$SCREEN_ROOT/pilot/run"
```

Inspect `$SCREEN_ROOT/pilot/run/sites.tsv` and its boxes against the prepared
receptor. Confirm that the actual site count and ligand count match the plan.
Check the raw `docking/` directory for `-best.pdbqt` poses, and
`results/ultidock_results.db` plus any `*-docking-results.csv` exports:

```bash
find "$SCREEN_ROOT/pilot/run/docking" -name '*-best.pdbqt' | wc -l
ls "$SCREEN_ROOT/pilot/run/results"
du -sh "$SCREEN_ROOT/pilot/run" "$SCREEN_ROOT/pilot/receptors"
```

Review failures and missing poses instead of silently dropping them. A filtered
CSV may be empty even when the database contains finite scores; see
[Results and reports](../user-guide/results-reports.md) for an unfiltered export.
Open representative best poses with the matching receptor. A site or preparation
problem found here should be fixed before a full library run.

## 5. Measure and tune a pilot

The repository does not provide a hardware-neutral ligands-per-hour guarantee.
Use the pilot's **completed ligand–site jobs per wall hour**, observed memory and disk
growth on the machine that will run the real screen. A rough first estimate is
`full ligand–site pairs / measured completed pairs per hour`, plus setup,
staging and analysis time. Receptor/site/grid work is partly fixed, so this is
an estimate, not a deadline. Retry or failure rates also matter.
For example, 40 completed ligand–site jobs in one hour suggests roughly
100 hours for 4,000 pairs on the same machine with the same settings, before
allowing for changed chemistry, setup and failures.

| Change | Likely throughput or capacity effect | Check before accepting it |
| --- | --- | --- |
| More CPU cores and `ULTIDOCK_WORKERS` | More independent jobs can run at once until memory, I/O or CPU contention dominates | In CPU mode, keep workers × `--vina-cpu` within allocated cores; test the actual worker count. |
| Higher `--vina-cpu` | More threads per Vina job; fewer jobs may run concurrently | Compare completed pairs per hour, not one ligand's latency. |
| Higher `--vina-exhaustiveness` | More Vina search effort, usually longer per pair | Hold it constant across a ranked screen and compare sampling quality on controls. |
| GPU backend and `GPU_SLOTS_PER_DEV` | GPU search can improve throughput; extra slots increase simultaneous GPU jobs | Verify a complete GPU run and watch device memory; CPU and GPU scores are not interchangeable. |
| More sites or larger boxes | More ligand–site jobs and larger grids, with more memory/disk use | Inspect biological coverage and measure a pilot with the intended site count. |
| Faster local scratch and enough free space | Less waiting on large map, pose and database I/O | Measure `du -sh`, `df -h` and job logs; retain outputs before scratch cleanup. |

The CPU worker default is based on the **host** CPU count divided by
`VINA_CPU`, which can exceed a scheduler allocation. For example, an
eight-core allocation with `--vina-cpu 2` can start with
`export ULTIDOCK_WORKERS=4`; bind the job to its allocated cores. GPU runs
default to two slots per detected device. Reduce slots if memory is constrained
and measure any increase; native `NUMWI` is a build setting, not a quick
per-run throughput switch. See [HPC / batch screening](../user-guide/hpc-batch-screening.md).

## 6. Run the full library

After accepting the pilot's site and settings, run the full prepared library
sequentially in its own output and receptor/grid directories. The command below
uses the same CPU protocol; use `--mode gpu` only if the pilot also used GPU
and its runtime passed a real run:

```bash
time ultidock cavity --autosites 1 --mode cpu --vina-cpu 2 \
  --vina-exhaustiveness 8 --vina-seed 42 --skip-wget \
  --macro-mol-dir "$SCREEN_ROOT/full/receptors" \
  --ligands-dir "$LIBRARY" \
  --output-dir "$SCREEN_ROOT/full/run"
```

Compare `sites.tsv` and the generated boxes with the pilot before treating
the two runs as one protocol. Save the actual generated `docking/config.py`
from the active workspace and the run's `run_config.yaml`; a later setup in
that workspace can replace `config.py`. Do not launch concurrent screens
against the same `ULTIDOCK_HOME`, because they share generated configuration.
For cluster arrays, isolate homes and output directories as shown in the
[HPC screening example](hpc-screening.md).

## 7. Audit and rank results

Check discovered ligand count, attempted/completed jobs, missing poses and failed
preparations against the input manifest. Planned ligand–site work is the valid
ligand count multiplied by the **actual** number of sites; failures can reduce
completed work, while analysis filters only change the exported rows. Keep
receptor, site, ligand state, engine and run/model identifiers together when
merging shards; names alone
can collide. Use the scored `docking_file` best pose rather than the
`ligand_file` input for visual inspection.

Inspect top poses, clashes and preparation diagnostics. When experimental labels
exist, evaluate enrichment separately; a favorable docking score by itself is
not proof of binding. See [Scoring and ranking](../scientific-background/scoring-ranking.md).
