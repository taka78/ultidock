# Site-Prediction Benchmarks

This directory is the new benchmark lane for comparing ligand-binding-site
prediction methods on DUD-E, COACH420, and HOLO4K.

Canonical method IDs are:

- `cav-emps`: Ultidock's CaV-EMPS method, Cavity detection via Electrostatic Map Pocket Scoring.
- `fpocket`
- `p2rank`

The older `benchmarks/cavity_recovery_benchmark.py` stays in place for the
DUD-E docking-box recovery work. This suite is narrower and more paper-facing:
every method produces predicted pocket centers, and the evaluator compares those
centers against ligand atoms using the same style of DCA metrics used in
COACH420/HOLO4K reports.

## Why This Exists

DUD-E was useful for proving Ultidock can recover docking boxes without using
the co-crystallized ligand. COACH420 and HOLO4K are better for a paper claim
about predictive binding-site performance because they are established
site-prediction datasets, and published P2Rank comparisons report Top-n and
Top-(n+2) DCA success rates on them.

## Dataset Contract

Normalized datasets use this shape inside the selected benchmark output
directory:

```text
<output-dir>/normalized/<dataset>/<target_id>/
  receptor_input.pdb
  labels.json
  metadata.json
```

`receptor_input.pdb` is the protein structure that each method receives. It
must not expose the answer to methods that would treat ligand HETATM records as
ordinary input.

`labels.json` stores the reference binding sites:

```json
{
  "binding_sites": [
    {
      "site_id": "L1",
      "ligand_id": "ATP:A:401",
      "center": [1.0, 2.0, 3.0],
      "ligand_atoms": [[1.0, 2.0, 3.0], [1.5, 2.0, 3.2]],
      "source": "original pdb ligand atoms"
    }
  ]
}
```

The ligand center is retained for diagnostics, but the primary COACH420/HOLO4K
metric is DCA: distance from a predicted pocket center to the closest atom of
the reference ligand.

Reference ligands follow the P2Rank parameters distributed with these datasets:

- select one deterministic alternate-location conformer;
- ignore `HOH`, `DOD`, `WAT`, `NAG`, `MAN`, `UNK`, `GLC`, `ABA`, `MPD`, `GOL`,
  `SO4`, and `PO4` HET groups;
- merge HET groups connected within 1.7 A into one ligand;
- require at least five heavy atoms;
- require a ligand atom within 4 A of a protein atom;
- require the ligand center within 5.5 A of a protein atom.

Structures with no relevant ligand after these rules remain in the manifest
with zero reference sites. They are reported as non-evaluable and do not enter
site-level success-rate denominators. `--labels-only` regenerates labels and
metadata while preserving existing normalized receptors.

## Prediction Contract

Each method writes one normalized `predictions.tsv`:

```text
target_id	method	rank	site_id	center_x	center_y	center_z	score	source
1abc	cav-emps	1	S1	0.0	1.0	2.0	0.42	centers.tsv
```

The DUD-E cavity-recovery and baseline runners also write this table now, so
their outputs can be converted into the new site-prediction evaluator without a
method-specific parser.

Site labels such as `S1` are just method-local names. Evaluation depends on
coordinates and rank/order, not on the label text. For `cav-emps`, the
publication rank should be produced by sorting the CaV-EMPS ranking score
descending and writing that as a separate rank column. The source `S*` labels
should not be renamed to pretend they are ranks.

## Primary Metrics

Primary paper metrics:

- DCA threshold: 4 A.
- Top-n: consider only the first `n` predicted pockets, where `n` is the number
  of true ligand sites in that target.
- Top-(n+2): same, but with two extra predictions.
- Success rate: fraction of reference ligand sites recovered within the DCA
  threshold.

Secondary diagnostics:

- DCA at 5 A and 10 A.
- Per-target failures.
- Runtime and number of emitted predictions.
- DUD-E legacy crystal-center distance, kept only for continuity with old runs.

## External Dataset Source

For COACH420 and HOLO4K, Ultidock uses the `rdk/p2rank-datasets` repository
because it is the dataset source used by P2Rank/PrankWeb comparisons. The
benchmark command downloads that repository automatically when it is not already
available locally.

Reference URLs:

- `https://github.com/rdk/p2rank-datasets`
- `https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0285-8/tables/3`

## Command Flow

The normal researcher-facing command is the automatic benchmark runner:

```bash
ultidock benchmark site-prediction \
  --datasets coach420,holo4k \
  --methods cav-emps,fpocket,p2rank \
  --jobs 4 \
  --output-dir benchmarks/results/site_prediction/coach_holo_$(date +%Y%m%d_%H%M)
```

That single command downloads `rdk/p2rank-datasets` if needed, normalizes the
selected datasets, runs each selected method, evaluates DCA Top-n/Top-(n+2), and
writes HTML/Markdown reports.

For `cav-emps`, the runner uses Ultidock's bundled AutoGrid source and will
compile/check `docking/AUTODOCK_GPU_DIR/autogrid/autogrid4` automatically when
it is missing. `fpocket` and `p2rank` still require their external executables
(`fpocket` and `prank`) to be installed or passed with `--fpocket-cmd` /
`--p2rank-cmd`.

CaV-EMPS AutoGrid maps are treated as debug artifacts. By default they are
created under a temporary work root and deleted after `centers.tsv` and
`predictions.tsv` are written. Add `--keep-artifacts` to preserve per-target
`maps/` folders, or pass `--work-root /path/to/scratch` to put temporary work on
a larger disk.

For smoke testing, keep it tiny:

```bash
ultidock benchmark site-prediction \
  --datasets coach420 \
  --methods cav-emps \
  --max-targets 3 \
  --jobs 1 \
  --dry-run
```

Advanced commands are still available if a step needs to be rerun manually:

```bash
ultidock benchmark site-import ...
ultidock benchmark site-run ...
ultidock benchmark site-evaluate ...
ultidock benchmark site-report ...
```
