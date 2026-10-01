# Baseline Binding-Site Finders

This directory contains receptor-only baseline benchmarks for external binding-site
finders. Each method lives in its own folder and writes outputs into a matching
method folder under `benchmarks/results/baseline/`.

The evaluation methodology is intentionally tied to Ultidock's DUD-E
cavity-recovery benchmark: DUD-E receptor input only, no crystal ligand given to
the site finder, then post-hoc comparison of predicted site centers against the
held-out crystal ligand.

Outputs include two metric families:

- Legacy DUD-E docking-box recovery: distance from predicted site center to the
  crystal-ligand geometric center.
- Publication-style site prediction: DCC, measured from predicted site center to
  the closest crystal-ligand atom, with Rank-1, Rank-3, and any-site success
  fields.

For fpocket and P2Rank, `rank_by_method` follows the method's reported
score/order. For CaV-EMPS, the matching field is `rank_by_fitness`, derived from
the ranking score in `centers.tsv`. Site labels such as `S1` remain names, not
ranks.

By default, baselines receive the same receptor model as Ultidock: the DUD-E
receptor is prepared through molguard/Meeko into a canonical PDBQT, then converted
back to PDB coordinate records because P2Rank and fpocket consume PDB input. Use
`--baseline-receptor-model original` only when you intentionally want a
third-party-tool-on-raw-DUD-E comparison.

## P2Rank

```bash
OUT="benchmarks/results/baseline/p2rank_$(date +%Y%m%d_%H%M)"

ultidock benchmark baseline-p2rank \
  --dataset-root benchmarks/datasets \
  --targets all \
  --autosites 6 \
  --jobs 4 \
  --force \
  --output-dir "$OUT"
```

If `prank` is not on `PATH`, provide the executable:

```bash
ultidock benchmark baseline-p2rank \
  --p2rank-cmd "/path/to/prank predict" \
  --dataset-root benchmarks/datasets \
  --targets all
```

## fpocket

```bash
OUT="benchmarks/results/baseline/fpocket_$(date +%Y%m%d_%H%M)"

ultidock benchmark baseline-fpocket \
  --dataset-root benchmarks/datasets \
  --targets all \
  --autosites 6 \
  --jobs 4 \
  --force \
  --output-dir "$OUT"
```

Both runners produce:

- `summary.csv`
- `sites.csv`
- `predictions.tsv`
- `summary.json`
- `run_metadata.json`
- per-target `centers.tsv`
- per-target raw method output under `<target>/raw/<method>/`

Use `--dry-run` to confirm target selection before running a full DUD-E pass.
