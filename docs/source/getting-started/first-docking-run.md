# First Docking Run

Start with the bundled D2 example. It has one receptor and three prepared ligands,
so you can follow the whole calculation without downloading a screening library.

```bash
ultidock example list
ultidock example run d2-antipsychotics --dry-run
ultidock example run d2-antipsychotics --mode auto
```

The preview lists the inputs without creating a workspace. The real run stages
`6CM4-edited.pdbqt`, `haloperidol.pdbqt`, `escitalopram-e.pdbqt`, and
`morphine-e.pdbqt` in `examples/d2-antipsychotics/workspace/<timestamp>/`.
The runner skips ligand downloads and uses only those three ligands.

Watch four phases: receptor checks/preparation, binding-site discovery, grid
preparation, and ligand docking. If receptor recovery emits warnings, save its
`.prep.json` report with the results. A successfully parsed file does not establish
correct protonation or bond orders.

Within the printed workspace:

| Directory | What to inspect |
| --- | --- |
| `MACRO_MOL_DIR` | Prepared receptor, site folders, grid inputs and logs. |
| `LIGANDS_DIR` | The staged ligand inputs. |
| `DOCKING_DIR` | Engine outputs and best-pose PDBQT files. |
| `RESULTS_DIR` | SQLite database, site metadata, and filtered CSV exports. |

In a CSV row, `docking_file` is the saved best-pose output and `ligand_file` is the
input. Open the pose with the corresponding receptor, and verify that the reported
site/model matches what you are viewing. An empty filtered CSV is possible:
thresholds can exclude all rows even when docking completed successfully.

The `--mode auto` option uses a detected compatible GPU for AutoDock-GPU and
CPU Vina when no GPU is visible. Without this option, the D2 example defaults
to CPU mode. With several NVIDIA GPUs, jobs are distributed across the detected
devices. Use `--mode gpu` if GPU execution is required. See [Results & reports](../user-guide/results-reports.md) for
unfiltered analysis and [the D2 tutorial](../tutorials/small-molecule-docking.md)
for interpreting this comparison set.
