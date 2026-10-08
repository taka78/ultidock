# Results and reports

Keep receptor, ligand, binding site, scoring engine and model/run identity with
its output coordinates. The pipeline writes scored models to SQLite as docking
finishes, then creates filtered analysis exports.

## Find the outputs

| Location | Contents |
| --- | --- |
| `DOCKING_DIR` | Vina output PDBQTs or AutoDock-GPU XML/DLG files, pose provenance, run manifests and failure reports |
| `RESULTS_DIR/ultidock_results.db` | Successfully parsed scored models across the run's receptors, ligands and sites |
| `RESULTS_DIR/*-docking-results.csv` | Filtered analysis tables |
| `MACRO_MOL_DIR` | Prepared receptors, site definitions, grids and grid logs |

The high-level `known-site`, `cavity`, `blind`, `fpocket` and `p2rank` commands
create timestamped folders under `docking/RESULTS_DIR/`; `--output-dir` chooses
another location. A cavity run's database is
`docking/RESULTS_DIR/cavity_<timestamp>/results/ultidock_results.db` in the active
workspace. Its raw engine outputs are in the run's `docking/` subdirectory.
Example runners instead print a workspace with the uppercase directory names
above. Raw `ultidock run` uses its generated configuration paths.

One invocation stores successful records for all its receptor–ligand–site
combinations in one database. Separate high-level invocations create separate
databases; there is no automatic cross-run SQLite merge. Pocket-method runs
also retain their staged receptor, sites and predictor outputs.

## Match a row to its coordinates

| Field | Meaning |
| --- | --- |
| `ligand_name` | Combined docking-case name, including the receptor stem |
| `binding_site` | Site identifier; may be blank in older or unknown-site records |
| `model` | Vina MODEL identifier or AutoDock-GPU XML/DLG run identifier |
| `docking_file` | Vina output PDBQT container or AutoDock-GPU DLG with coordinates |
| Score and RMSD fields | The parsed model's values; unavailable GPU RMSD values remain missing |

The CSV derives `model` from the stored `ligand_name`; it is not a separate
model column in SQLite. This branch stores all successfully parsed scored models. Analysis applies its
score, RMSD and model filters; it does not automatically reduce each case to
one best-pose row. There is no `ligand_file` or `is_best_pose` database column
in this branch. Retain the staged inputs and provenance records to trace input
identity. Derived CSV identifiers are not a replacement for those records.

For Vina, select the row's MODEL from the output PDBQT. For AutoDock-GPU, inspect
the same run in its DLG and keep the accompanying score XML. The engine can also
write a `-best.pdbqt`, but that file's selection criterion can differ from the
lowest XML binding energy. It does not replace a specific row's run identity.
See [Scoring and ranking](../scientific-background/scoring-ranking.md).

## Analysis filters and reports

Headers are retained when filters exclude every pose. To export all stored
finite scores from a source checkout with a generated docking configuration:

```bash
python docking/analyse_docking_results.py --db /absolute/results/ultidock_results.db \
  --affinity inf --rmsd-lb inf --rmsd-ub inf --min-model 0 \
  --out /absolute/results/all-poses.csv
```

The minimum model filter is exclusive. Missing GPU RMSD bounds are accepted by
the analysis; they are not zero-distance measurements. Excel output requires
an `.xlsx`/`.xls` output path and the appropriate writer package.

`ultidock report /absolute/run-directory` builds Markdown/HTML and visualization
helpers from recognized files such as `sites.tsv`, `predictions.tsv`,
`top_hits.csv` or `scores.csv`. It does not automatically turn every SQLite
schema or timestamped CSV into a populated hit table. Inspect engine exports
and SQLite directly if a generic report section is empty.

When sites are available, helpers can include `cavity_centers.pdb`,
`cavemps_sites.pml`, `site_boxes.pml`, `top_poses.pml`, `cavemps_sites.cxc` and
`cavity_centers.bild` for PyMOL/ChimeraX inspection. Archive configuration,
source structures, prepared inputs, site definitions, seeds, tool versions,
logs and raw outputs with derived figures or hit lists.

## Partial screens and MD handoff

Each pipeline invocation writes a `docking-run-*.json` manifest of its
successful outputs, including hashes and identities, plus a sibling
`docking-run-*.failures.csv` with stage, receptor, ligand, site and error.
A failed preparation or docking case is reported while other cases continue.
Partial screens proceed into analysis and optional MD. If no docking case
succeeds, the run returns failure after reporting and skips those stages.
Shared setup/backend failures can still stop the pipeline. Generated reports
include a Screening outcomes section with status, successful-output count and
a failure-table preview; retain the full failure CSV for complete accounting.

The `.md.json` provenance sidecars beside pose outputs bind them to receptor
and ligand inputs. An integrated MD handoff uses the current run manifest so
old results in a reused directory cannot enter its selection. The sibling
`docking-run-*.md.json` records the MD handoff status and job path.

MD selects raw scored poses independently of the CSV filters, for the protocol's
one receptor and engine. Its systems, checkpoints, trajectories and energy files
live under `md-simulation/workspace/`, not in the docking SQLite database.
See [Docking to MD](molecular-dynamics/docking-to-md.md) and
[MD execution and results](molecular-dynamics/running-and-results.md).
