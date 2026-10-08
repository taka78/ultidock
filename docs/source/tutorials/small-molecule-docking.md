# Small-molecule docking

Use the three-ligand 4COF example to follow receptor checks, site generation,
docking and scored-pose inspection. Complete [Installation](../getting-started/installation.md)
first. GPU execution needs its driver/runtime and AutoDock-GPU; automatic mode
uses CPU Vina when no GPU is available.

## 1. Inspect the supplied inputs

```bash
ultidock example run gabaa-benzos --dry-run
```

The runner stages `4COF_edited.pdbqt`, `aspirine-e.pdbqt`, `ibuprofen-e.pdbqt`
and `morphine-e.pdbqt` in a new workspace, then prints the planned command.
These are the branch's actual inputs despite the legacy `gabaa-benzos` name.
No other files from your usual ligand directory are added.

## 2. Run docking

```bash
ultidock example run gabaa-benzos
```

The default automatic mode prefers a compatible GPU and supports multiple
visible NVIDIA devices. Add `--mode gpu` to require AutoDock-GPU. Ultidock
checks the receptor, proposes CaV-EMPS sites, builds grids, schedules all three
ligands and writes SQLite records and analysis exports automatically.

## 3. Inspect a scored model

Open the CSV in the printed workspace's `RESULTS_DIR/`. Follow `docking_file`
to the Vina PDBQT container or AutoDock-GPU DLG in `DOCKING_DIR/`, and inspect
the row's particular model/run with the matching prepared receptor. Check the
site, coordinates, contacts and clashes. The ligand in `LIGANDS_DIR/` is an
input, rather than the engine's sampled output.

An empty CSV may mean that all poses were filtered out. Use the unfiltered
export in [Results and reports](../user-guide/results-reports.md) to inspect
all stored models. Scores alone do not label these comparison molecules as
validated binders or decoys.

## 4. Compare another prediction method

For a bundled runner that supports multiple site finders, use SERT:

```bash
ultidock example run sert-escitalopram --mode gpu
ultidock example run sert-escitalopram fpocket --mode gpu
```

Install fpocket first, then compare the new workspaces' boxes and pose outputs.
Keep the receptor, ligand and docking engine fixed when assessing the change
of site method. See [Bundled examples](examples.md) for available runners and
[Docking to MD](../user-guide/molecular-dynamics/docking-to-md.md) for subsequent
simulation with a reviewed protocol.
