# Results & reports

Keep input identity, receptor identity, site, engine run/model and pose file
together. A ligand input is not a docking result, even if it shares the molecule's
name. Example runs print an isolated workspace containing:

- `DOCKING_DIR`: raw engine outputs and best-pose PDBQT files.
- `RESULTS_DIR/ultidock_results.db`: scored runs and associated metadata.
- `RESULTS_DIR/*-docking-results.csv`: filtered analysis exports.
- `MACRO_MOL_DIR`: prepared receptors, grids, sites and preparation reports.

New analysis exports select the saved best pose per ligand/receptor/site/run
invocation. `docking_file` points to that PDBQT; `ligand_file` records the prepared
input. All parsed runs remain in SQLite. Legacy rows lacking best-pose metadata
remain exportable, so mixed old/new databases need care.

The default analysis filters can produce an empty CSV. To inspect all finite
scores from a database in a source checkout:

```bash
python docking/analyse_docking_results.py --db /absolute/results/ultidock_results.db \
  --affinity inf --rmsd-lb inf --rmsd-ub inf --out /absolute/results/best-poses.csv
```

The analysis script imports a generated docking configuration. Run setup first
or use the configuration from the completed workflow. `model` is the engine's
run/model identifier, not a site rank. Missing GPU RMSD bounds are not zeros.

`ultidock report /absolute/run-directory` writes Markdown/HTML and visualization
helpers from recognized files such as `sites.tsv`, `predictions.tsv`,
`top_hits.csv` or `scores.csv`. The generic report does not automatically convert
every SQLite database or timestamped CSV into a populated top-hits table. Inspect
the engine CSV/database directly if a report section is empty.

Open a best PDBQT and its matching prepared receptor in a molecular viewer.
Inspect placement, clashes, interactions and preparation warnings before choosing
candidates for further work. See [Scoring / ranking](../scientific-background/scoring-ranking.md).
