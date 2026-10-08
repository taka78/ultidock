# Run the bundled examples

Examples let you run Ultidock before supplying your own molecules. They use
bundled receptor and prepared ligand PDBQT files; dataset examples also include
`dataset.json` metadata describing their inputs. This is separate from
`docking/ligands.wget`, the optional download list for your own runs. Example
scripts pass `--skip-wget` so that download list cannot add molecules to the
example.

## Try one complete run

From an installed checkout with the required tools available:

```bash
ultidock example list
ultidock example run d2-antipsychotics --dry-run
ultidock example run d2-antipsychotics --mode auto
```

`example list` shows the runnable examples. `--dry-run` previews the selected
files without docking. The real command creates a new
`examples/d2-antipsychotics/workspace/<timestamp>/` folder, copies the declared
inputs, prepares the receptor, finds sites, builds grids and docks the ligands.
It prints the workspace path. Existing example runs remain available for
inspection; installed native tools are reused.

Open `RESULTS_DIR/` inside that workspace for the SQLite score database and
CSV tables. The corresponding 3D output poses are in `DOCKING_DIR/`. Follow
[First Docking Run](../getting-started/first-docking-run.md) to inspect a pose
and [Results and reports](../user-guide/results-reports.md) for the output
files.

## Choose another example

| Example | What it demonstrates | Command |
| --- | --- | --- |
| `quickstart` | Interactive guided D2 run | `ultidock example run quickstart` |
| `d2-antipsychotics` | Dopamine D2 receptor with three comparison ligands | `ultidock example run d2-antipsychotics --mode auto` |
| `sert-escitalopram` | Serotonin transporter with one ligand | `ultidock example run sert-escitalopram --mode auto` |
| `gabaa-benzos` | GABA-A receptor and ligand workflow | `ultidock example run gabaa-benzos --mode auto` |
| `gabaa-8dd2-cav-emps` | Site prediction against withheld known sites | `ultidock example run gabaa-8dd2-cav-emps --dry-run` |

The D2 example stages `6CM4-edited.pdbqt`, haloperidol, escitalopram and
morphine. Its [dataset notes](https://github.com/taka78/ultidock/tree/v1.1.2/examples/d2-antipsychotics)
give input provenance; the comparison molecules are not validated decoys.
The GABA 8DD2 study withholds five GABA/zolpidem sites and compares controlled
CaV-EMPS settings with fpocket and P2Rank. Read its example README before
interpreting the results.

Bundled examples have different defaults: D2 and the interactive teacher
start in CPU mode, while SERT and GABA-A require a GPU by default. The
`--mode auto` commands above select a visible compatible GPU for AutoDock-GPU
or fall back to CPU Vina; the teacher accepts only `cpu` or `gpu`.
On a CUDA machine with multiple NVIDIA GPUs, docking jobs are spread across
detected devices. Pass `--mode gpu` when GPU execution is required.

The `quickstart` teacher explains each stage before launching a real D2 run.
Press **Enter** to advance, **b** to go back and **q** to exit. Preview all
lessons with `ultidock example run quickstart --dry-run`. Use `--yes` only for
explicit unattended execution. See [Quick Start](../getting-started/quick-start.md).

## Try another site finder

D2 and SERT use Ultidock's included CaV-EMPS finder by default. Pass `fpocket`
or `p2rank` after the example name to compare a different prediction method:

```bash
ultidock example run sert-escitalopram fpocket --mode auto
ultidock example run sert-escitalopram p2rank --mode auto
```

These methods need their optional local tools; see
[Installation](../getting-started/installation.md). P2Rank 2.5 needs Java
17–23; the documented Ubuntu choice is `openjdk-21-jre-headless`. If Java 25
is selected, it may fail with `Unsupported class file major version 69`.
The bundled `external/bin/prank` launcher uses
`/usr/lib/jvm/java-21-openjdk-amd64` when `JAVA_HOME` is unset; set
`JAVA_HOME` to Java 21 if it points to an incompatible installation.

The predictor writes new site boxes before docking. Each site needs its own
AutoGrid map set; `[autogrid] prepared N site grids` means grid preparation
finished. Docking work grows with the number of ligands times the number of
sites. Keep the chosen method, tool version and raw prediction output with
comparisons. The example scripts can also serve as automation templates.
