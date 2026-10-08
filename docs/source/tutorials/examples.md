# Run the bundled examples

Examples supply inputs and runnable scripts so you can try Ultidock before
staging your own molecules. Discover the examples available in this branch:

```bash
ultidock example list
ultidock example run sert-escitalopram --dry-run
ultidock example run sert-escitalopram
```

SERT and the 4COF docking example create a fresh timestamped workspace, copy
the declared receptor and ligand PDBQTs, and invoke the full pipeline. Their
scripts disable `docking/ligands.wget`, so a download list cannot add unrelated
ligands. Their dry runs stage the files and print the command without docking.
Use the printed workspace to find `RESULTS_DIR/ultidock_results.db`, filtered
CSV tables and `DOCKING_DIR/` engine outputs.

## Choose an example

| Example | Included inputs and purpose |
| --- | --- |
| `sert-escitalopram` | `5i6x_edited.pdbqt` and `escitalopram-e.pdbqt`; native docking and alternative site finders |
| `gabaa-benzos` | `4COF_edited.pdbqt` with `aspirine-e.pdbqt`, `ibuprofen-e.pdbqt` and `morphine-e.pdbqt`; three-ligand docking workflow |
| `gabaa-8dd2-cav-emps` | Receptor-only recovery of withheld GABA/zolpidem sites; a structural site-prediction study |
| `quickstart` | Illustrative report/artifact generation, with no native docking or MD |

The legacy name `gabaa-benzos` does not describe a benzodiazepine panel in this
branch; use the filenames and dataset notes to identify its actual inputs.
The `examples/d2-antipsychotics/` directory contains a receptor and notes, but
has no runnable example or bundled ligand panel here.

Read each example's README and `dataset.json` when provided for provenance and
limitations. Comparison molecules are not automatically validated actives or
decoys. The [small-molecule walkthrough](small-molecule-docking.md) uses the
three-ligand 4COF example.

## Use the GPU backend

SERT and 4COF default to `--mode auto`: detected GPU execution is preferred,
with Vina as CPU fallback. CUDA docking uses multiple visible NVIDIA GPUs.
Require GPU execution when comparing AutoDock-GPU results:

```bash
ultidock example run sert-escitalopram --mode gpu
ultidock example run gabaa-benzos --mode gpu
```

The report demonstration has no hardware mode. See [Quick Start](../getting-started/quick-start.md)
for its outputs and [First Docking Run](../getting-started/first-docking-run.md)
for inspecting native results.

## Compare site finders

SERT accepts an optional positional site method. CaV-EMPS is included; fpocket
and P2Rank require their [optional installations](../getting-started/installation.md).

```bash
ultidock example run sert-escitalopram fpocket --mode gpu
ultidock example run sert-escitalopram p2rank --mode gpu
```

Each prediction generates site boxes and one AutoGrid map set per site before
docking. Keep settings, predictor versions and raw predictions with comparisons.
P2Rank 2.5 needs Java 17–23; Java 21 is the documented Ubuntu choice. See
[P2Rank](../user-guide/binding-site-discovery/p2rank.md) for Java troubleshooting.

## Continue into MD

The SERT and 4COF runners accept `--md-config`, `--md-through` and `--md-work-dir`.
However, their PDBQT docking inputs do not supply a reviewed protein PDB,
complete ligand SDF/atom map or membrane seed. Complete those inputs using the
[MD guide](../user-guide/molecular-dynamics/index.md) before adding a protocol.
The default continuation ends at NPT; production requires an equilibration review.
