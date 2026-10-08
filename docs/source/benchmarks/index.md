# Benchmarking

Benchmark scripts are in `benchmarks/`. Downloaded DUD-E datasets and generated
results are ignored by Git and should remain local. Use
[System requirements](../getting-started/requirements.md) to provision native
tools before benchmarking. A site-recovery benchmark tests localization;
active/decoy docking evaluates ranking and is much more expensive.

## Download DUD-E targets

The download command fetches receptor, co-crystal ligand, active and decoy
files for selected targets:

```bash
ultidock benchmark download-dude \
  --targets ace,bace1,braf,cdk2,cxcr4,drd3,egfr,esr1,gcr,hdac2,hivpr,pde5a,pparg,src,vgfr2 \
  --dataset-root benchmarks/datasets
```

The corresponding checkout script is `python3 benchmarks/download_dude.py`.
Keep dataset versions, target lists and preparation decisions with results.

## Receptor-only site recovery

```bash
ultidock benchmark cavity-recovery \
  --dataset-root benchmarks/datasets \
  --targets ace,bace1,braf,cdk2,cxcr4,drd3,egfr,esr1,gcr,hdac2,hivpr,pde5a,pparg,src,vgfr2 \
  --autosites 6 --site-policy receptor_search --jobs 4 --force \
  --output-dir benchmarks/results/cavity_recovery
```

The benchmark predicts sites without the crystal ligand and evaluates them
afterward. It records both center-to-centroid distance and DCC distance to
the nearest ligand atom. These measure different localization questions.
The output includes `summary.csv` (target-level closest-site, DCC and Top-k
flags), `sites.csv` (per-site distances and ranking metadata),
`predictions.tsv` (normalized CaV-EMPS predictions), and each target's
`centers.tsv` boxes. The underlying script is
`benchmarks/cavity_recovery_benchmark.py`; `--jobs` parallelizes targets.

Large AutoGrid maps and scratch files are temporary by default. Use
`--keep-artifacts` for a debugging target and `--work-root` to place scratch
on a larger disk. Do not interpret a site ID such as S1 as a quality rank;
report all proposals and the declared Top-k budget.

## COACH420 and HOLO4K site prediction

```bash
ultidock benchmark site-prediction \
  --datasets coach420,holo4k \
  --methods cav-emps,fpocket,p2rank \
  --jobs 4 \
  --output-dir benchmarks/results/site_prediction/coach_holo_$(date +%Y%m%d_%H%M)
```

The command downloads `rdk/p2rank-datasets` when needed, normalizes receptors
and ligand labels, runs each chosen predictor, evaluates DCC Top-n and
Top-(n+2), and writes Markdown/HTML reports. Prediction and evaluation are
separate stages: reference ligand geometry is not fed to receptor-only
prediction. The CaV-EMPS arm uses compact scratch by default;
`--keep-artifacts` preserves large maps for debugging and
`--work-root /path/to/scratch` moves scratch to a larger drive.
See [Binding-site prediction](../scientific-background/binding-site-prediction.md)
for evaluation design.

## Active/decoy docking and enrichment

`benchmarks/dude_docking_benchmark.py` and
`benchmarks/run_full_dude_benchmark.py` run full active/decoy docking arms.
They collect ROC-AUC, EF1, EF5, BEDROC and logAUC. These are final validation
runs with substantially greater compute and storage cost than receptor-only
site recovery. Keep failed targets and missing ligands in the accounting and
do not compare Vina and AutoDock-GPU raw scores as one scale.

The [results guide](../user-guide/results-reports.md) explains saved poses and
exports, while the [screening walkthrough](../tutorials/virtual-screening.md)
shows how to measure a small pilot before a large job.
