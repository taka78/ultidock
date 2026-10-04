# GABAA 4COF Mixed-Ligand Docking Example

This example stages the prepared 4COF receptor and three bundled ligands:
`aspirine-e`, `ibuprofen-e` and `morphine-e`. The directory retains its legacy
`gabaa-benzos` name; these files are not a benzodiazepine panel. Set
up Python, AutoGrid, and a docking backend using the [setup guide](../../SETUP.md)
before starting a full run.

From the repository root:

```bash
ultidock example run gabaa-benzos
# On a machine without a GPU runtime:
ultidock example run gabaa-benzos --mode cpu
```

For a source-checkout installation, replace `ultidock` with
`/usr/bin/python3 -m cli.ultidock`. The runner creates a new directory under
`workspace/<timestamp>/`, copies only the bundled receptor and ligands, skips
the general ligand download manifest, and prints the directory path so you can
inspect the results. You can also run `python3 example-run.py` from this folder.

The default `--mode auto` permits CPU fallback. `--mode gpu` explicitly requires
a GPU. `--dry-run` stages the inputs and prints the pipeline command;
`--output-dir ./runs/gabaa` selects a new workspace relative to your current
directory. Existing output directories are rejected to preserve earlier runs.

Receptor numeric-format warnings are corrected on the staged copy by the
normal preparation step. Individual screening failures are logged in
`DOCKING_DIR/docking-run-*.failures.csv` while successful cases continue.

`--md-config /path/to/protocol.json` enables the shared MD handoff, with the
same `--md-through`, `--md-work-dir` and MD tool overrides as `ultidock run`.
The bundled inputs are docking-only: membrane MD also requires a reviewed
protein PDB, ligand SDFs/atom maps and compatible bilayer/orientation inputs.
Use `receptor_id=4COF_edited` and the ligand keys listed above; see the
[MD guide](../../md-simulation/README.md).

Structure source: [PDB 4COF](https://www.rcsb.org/structure/4COF).
