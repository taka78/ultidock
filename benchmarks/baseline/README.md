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
receptor is prepared through molguard with Meeko or Open Babel into a canonical PDBQT, then converted
back to PDB coordinate records because P2Rank and fpocket consume PDB input. Use
`--baseline-receptor-model original` only when you intentionally want a
third-party-tool-on-raw-DUD-E comparison.

Install the Python and native prerequisites in [SETUP.md](../../SETUP.md), then
run `bash scripts/install_pocket_tools.sh fpocket` and
`bash scripts/install_pocket_tools.sh p2rank` from the repository root. The
baseline commands do not install these tools automatically, so pass the local
binary paths shown below. P2Rank 2.5 needs Java 17–23; the Ubuntu package list
installs Java 21.

## P2Rank

```bash
OUT="benchmarks/results/baseline/p2rank_$(date +%Y%m%d_%H%M)"

ultidock benchmark baseline-p2rank \
  --dataset-root benchmarks/datasets \
  --targets all \
  --autosites 6 \
  --jobs 4 \
  --force \
  --p2rank-cmd "$PWD/external/bin/prank predict" \
  --output-dir "$OUT"
```

If using another P2Rank installation, provide its executable:

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
  --fpocket-cmd "$PWD/external/fpocket/bin/fpocket" \
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
