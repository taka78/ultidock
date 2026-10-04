# Bundled examples

List runnable entries with `ultidock example list`. Each docking example uses
an isolated workspace so repeated runs cannot mix inputs or scores. Tools are
shared with the active Ultidock workspace.

| Example | What it runs | Requirements beyond the Python package |
| --- | --- | --- |
| [quickstart](quickstart/readme.md) | Synthetic artifact and report demonstration; no docking or MD | None |
| [sert-escitalopram](sert-escitalopram/readme.md) | One SERT receptor and escitalopram, with CaV-EMPS, P2Rank or fpocket sites | AutoGrid plus Vina or AutoDock-GPU; selected pocket tool and its build/runtime dependencies |
| [gabaa-benzos](gabaa-benzos/readme.md) | 4COF docking with aspirin, ibuprofen and morphine inputs | AutoGrid plus Vina or AutoDock-GPU |
| [gabaa-8dd2-cav-emps](gabaa-8dd2-cav-emps/readme.md) | Preregistered site-recovery study; no ligand docking or MD | Recorded 8DD2 structure (local or downloaded), AutoGrid, receptor preparation tools, and the selected pocket predictors |
| [d2-antipsychotics](d2-antipsychotics/readme.md) | Receptor data only; no runnable entry point | Add ligand inputs and use the main pipeline |

Start with `ultidock example run quickstart`. The SERT and 4COF examples default
to `--mode auto`, permitting CPU fallback. Their `--dry-run` option stages
inputs and prints the command; `--output-dir PATH` chooses a new workspace.
AutoGrid can be built during setup when its compiler prerequisites are present.
Run `ultidock doctor` to see which tools are available before a full run.

The membrane-receptor examples are substantial CPU jobs. CaV-EMPS first builds
whole-receptor AutoGrid maps, then selects sites and builds their docking grids.
Map generation can take several minutes or longer before docking starts. The
pipeline prints the absolute `grid.glg` path; use `tail -f /path/to/grid.glg`
in another terminal to follow AutoGrid's progress. A short timeout does not
establish whether the full example works.

SERT and 4COF accept `--md-config PROTOCOL` and the other MD continuation flags.
**None of the examples bundles a complete, reviewed MD system.** Their PDBQT
files alone do not provide ligand bond orders/atom maps or membrane parameters.
Complete the [MD inputs](../md-simulation/README.md) before requesting simulation.

Individual docking failures are skipped and recorded in the run manifest and
its `.failures.csv` report. Successful cases continue through analysis and,
when requested and valid, MD. Missing shared docking infrastructure still needs
to be installed before screening can run.
