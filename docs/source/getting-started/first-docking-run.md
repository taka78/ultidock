# First Docking Run

The bundled SERT example supplies a serotonin-transporter receptor and one
prepared escitalopram ligand. Install [Ultidock and the native tools](installation.md)
first, then run:

```bash
ultidock example list
ultidock example run sert-escitalopram --dry-run
ultidock example run sert-escitalopram
```

The example stages `5i6x_edited.pdbqt` and `escitalopram-e.pdbqt` in a fresh
`examples/sert-escitalopram/workspace/<timestamp>/` folder. Its download stage
is disabled, so only the declared ligand is included. The dry run also stages
inputs and prints the planned pipeline command, but does not build grids or dock.
The real run prints its own workspace path.

The default automatic backend selects a detected compatible GPU for
AutoDock-GPU, with work distributed across visible NVIDIA GPUs on CUDA systems.
It falls back to CPU Vina if no GPU is detected. Require GPU execution with:

```bash
ultidock example run sert-escitalopram --mode gpu
```

Watch receptor checks, CaV-EMPS site discovery, AutoGrid preparation and docking.
Keep conversion warnings and engine logs. These checks validate formats and
execution; they do not establish correct biological preparation.

| Directory in the printed workspace | What to inspect |
| --- | --- |
| `MACRO_MOL_DIR` | Receptor, sites, grids and grid logs |
| `LIGANDS_DIR` | Prepared ligand input |
| `DOCKING_DIR` | Scored engine outputs, run manifest and failure report |
| `RESULTS_DIR` | SQLite docking database and filtered CSV exports |

A CSV row identifies a scored model/run and binding site. Its `docking_file`
points to the Vina output PDBQT container or AutoDock-GPU DLG. Inspect that
specific model with its matching receptor; the input ligand is not the docked
pose. Empty filtered exports can occur after successful docking. See
[Results and reports](../user-guide/results-reports.md) for interpretation.

For more examples and alternative site methods, use
[Run the bundled examples](../tutorials/examples.md). An optional
[docking-to-MD continuation](../user-guide/molecular-dynamics/docking-to-md.md)
requires reviewed chemical inputs and, for SERT, a complete membrane system.
The docking example alone does not supply those MD inputs.
