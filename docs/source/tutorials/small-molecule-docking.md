# Small-molecule docking

Run the bundled D2 receptor comparison to learn preparation, site generation and
best-pose inspection. Complete [Installation](../getting-started/installation.md)
first; GPU mode needs a compatible runtime, while CPU fallback needs Vina.
AutoGrid and receptor recovery dependencies are used in either case.

## 1. Inspect the inputs

```bash
ultidock example run d2-antipsychotics --dry-run
```

The manifest lists `6CM4-edited.pdbqt` as receptor and `haloperidol.pdbqt`,
`escitalopram-e.pdbqt`, and `morphine-e.pdbqt` as ligands. No other files from your
usual ligand directory are added. The prepared files are supplied example inputs;
retain their provenance when interpreting the run.

## 2. Dock the three ligands

```bash
ultidock example run d2-antipsychotics --mode auto
```

The runner creates a fresh timestamped workspace, copies the listed molecules,
prepares the receptor, proposes CaV-EMPS sites, builds grids and runs docking.
Review any hydrogen-recovery warnings. Original example molecules remain in the
example directory. For guided prompts, use `ultidock example run quickstart`.

## 3. Read the results

Follow the workspace printed in the terminal. Open the CSV in `RESULTS_DIR` and
follow `docking_file` to the selected best-pose PDBQT in `DOCKING_DIR`. Load that
pose with the prepared receptor from the same workspace. Check the site and model
columns; the ligand input in `LIGANDS_DIR` is not the docked pose.

If filtering yields no rows, inspect the database and relax analysis filters as
shown in [Results & reports](../user-guide/results-reports.md). Empty exports are
not evidence that no ligand can bind.

## 4. Change one variable

With fpocket installed, run:

```bash
ultidock example run d2-antipsychotics fpocket --mode auto
```

Compare search regions and poses across the two workspaces. Do not label morphine
or escitalopram as validated decoys based on this example. For provenance, see the
[dataset notes](https://github.com/taka78/ultidock/tree/v1.1.2/examples/d2-antipsychotics).
