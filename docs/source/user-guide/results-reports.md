# Results & reports

Keep input identity, receptor identity, site, engine run/model and pose file
together. A ligand input is not a docking result, even if it shares the molecule's
name. Example runs print an isolated workspace containing:

- `DOCKING_DIR`: raw engine outputs and best-pose PDBQT files.
- `RESULTS_DIR/ultidock_results.db`: scored runs and associated metadata.
- `RESULTS_DIR/*-docking-results.csv`: filtered analysis exports.
- `MACRO_MOL_DIR`: prepared receptors, grids, sites and preparation reports.

The high-level `known-site`, `cavity`, `blind`, `fpocket` and `p2rank`
commands create separate timestamped run folders under
`docking/RESULTS_DIR/` by default; `--output-dir` selects another folder.
For `ultidock cavity`, the run database is
`docking/RESULTS_DIR/cavity_<timestamp>/results/ultidock_results.db` in the
active workspace. One invocation writes successfully parsed records for all
its receptor–ligand–site combinations into that database. The
`ligand_name` includes the receptor stem, `binding_site` identifies the site,
and `ligand_file` and `docking_file` distinguish the input from the pose.
A separate invocation creates a separate run database; there is no automatic
cross-run SQLite merge.
Pocket-method folders also keep the staged receptor, predicted sites and raw
predictor output. Other modes keep receptor grids in the configured
`MACRO_MOL_DIR`. Older top-level `runs/` folders remain where they were.

New analysis exports select the saved best pose per ligand/receptor/site/run
invocation. `docking_file` points to that PDBQT; `ligand_file` records the prepared
input. All parsed runs remain in SQLite. Legacy rows lacking best-pose metadata
remain exportable, so mixed old/new databases need care.

The SQLite database receives incremental docking records, so it can be
inspected while a run is in progress. Analysis exports include
`binding_site` to distinguish identical ligand/model names from different
pockets; that field can be blank for unknown sites or older databases.
Exports retain headers even if filters exclude all poses. The best GPU pose
file is selected by AutoDock-GPU's total score and matched back to its DLG
coordinates for the corresponding score and run ID; it need not be the
lowest binding-energy XML row. Vina's best model is written separately.
Raw outputs and other parsed runs remain available in SQLite.
The analysis script writes CSV by default and supports Excel output when
given an `.xlsx` or `.xls` output path and the required writer is installed.

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

When site coordinates are present, the report generator can write
`report.md`, `report.html`, `cavity_centers.pdb`,
`cavemps_sites.pml`, `site_boxes.pml`, `top_poses.pml`,
`cavemps_sites.cxc` and `cavity_centers.bild`. The helper files
support PyMOL or ChimeraX inspection and retain the command, configuration,
software version and CaV-EMPS site ranking. Keep source structures,
prepared inputs, site definitions, seeds, engine versions, logs and raw poses
with any derived figures or hit list.

Open a best PDBQT and its matching prepared receptor in a molecular viewer.
Inspect placement, clashes, interactions and preparation warnings before choosing
candidates for further work. See [Scoring / ranking](../scientific-background/scoring-ranking.md).
