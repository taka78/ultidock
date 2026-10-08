# Quick Start: explore the outputs

This branch's quickstart creates a small demonstration report so you can learn
Ultidock's file layout without native docking or MD tools. Its site coordinates
and score rows are illustrative; they are not calculated docking results.

```bash
ultidock example run quickstart --dry-run
ultidock example run quickstart
```

The default output is `examples/quickstart/workspace/quickstart_run/` in the
active workspace. Open its `report.html`. Repeating the command writes to the
same default directory; use `--output-dir /absolute/path/to/demo` for another
location. The dry run only prints planned outputs.

| File | What you learn |
| --- | --- |
| `input/receptor.pdb`, `input/reference_ligand.mol2` | Small bundled structural fixtures |
| `run_config.yaml` | Command, workflow and version metadata |
| `sites.tsv`, `predictions.tsv` | Site coordinates and prediction table layout |
| `top_hits.csv` | Illustrative hit-table columns, with placeholder pose paths |
| `results.sqlite` | Demonstration metadata database, separate from the real docking database schema |
| `report.md`, `report.html` | Report presentation and links |
| Viewer helpers | Site-center and box inspection scripts |

The demonstration does not run a docking engine or create the placeholder
ligand poses. To test native docking and GPU execution, continue to
[First Docking Run](first-docking-run.md). To use your own inputs, follow
[Dock your own molecules](../user-guide/start-docking.md).
