# Dock your own molecules

Ultidock tests **ligands** (small molecules) against a **receptor** (the target
structure). The docking engines use **PDBQT**, a format that stores atom types,
charges and coordinates. A single `ultidock run` prepares the receptor, finds
search regions, builds energy maps called **grids**, docks the ligands and
analyzes the results.

## 1. Check your machine and install

Use Linux with Python 3.10 or later. CPU docking uses AutoDock Vina. GPU docking
needs a compatible NVIDIA CUDA or OpenCL device, driver and runtime. Check the
[system requirements](../getting-started/requirements.md), including the Ubuntu
native-package command, before selecting a backend. From a source checkout,
install the Python package in a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
ultidock doctor
```

`ultidock doctor` shows the active **workspace** (the folder holding workflow
scripts, inputs and results) and checks for installed tools. For a regular
package install or another Linux distribution, follow [Installation](../getting-started/installation.md).

## 2. Add a receptor

Put a `.pdb`, `.mol2` or `.pdbqt` receptor in the workspace's
`docking/MACRO_MOL_DIR/`. The pipeline prepares and checks it automatically.
In an editable checkout, the workspace is the repository by default; run these
commands from its root:

```bash
mkdir -p docking/MACRO_MOL_DIR
cp /path/to/receptor.pdb docking/MACRO_MOL_DIR/
```

For a regular install, put the workspace path printed by `ultidock doctor`
before `docking/`, for example `<workspace>/docking/MACRO_MOL_DIR/`. Choose
a receptor structure appropriate to the experiment;
automatic format preparation cannot decide its biological state for you.

## 3. Choose where ligands come from

Docking needs at least one ligand. You can select a chemically relevant slice
of a public library instead of assembling ligand files one by one. Ultidock
accepts the resulting AutoDock PDBQT download list directly:

| Source | What to provide | What `ultidock run` does |
| --- | --- | --- |
| ZINC 3D subset | A `wget` list for AutoDock `.pdbqt.gz` archives, saved as `docking/ligands.wget` | Downloads the selected archives, extracts and splits their ligands |
| Local files | Prepared `*.pdbqt` ligands or `.pdbqt.gz` archives in `docking/LIGANDS_DIR/` | Uses these inputs; add `--skip-wget` to avoid the bundled example download |

`wget` is a command-line file downloader. `ligands.wget` is a plain text list
of its commands. It lets a run fetch the same selected library again without
copying archives by hand. The repository's file contains one example download.
Each nonempty line must be a `wget` command, not a URL alone. For example:

```text
wget https://example.org/my-library.pdbqt.gz -O my-library.pdbqt.gz
```

Setup downloads listed `.pdbqt.gz` archives into `LIGANDS_DIR`. An archive
may contain many ligands; the pipeline unpacks it and uses Vina's `vina_split`
to make individual ligand PDBQT files.

### Start with a selected ZINC library

1. In the [ZINC20 3D tranche browser](https://zinc20.docking.org/tranches/home/)
   or [ZINC-22 tranche browser](https://cartblanche22.docking.org/), choose
   **3D** and select a manageable subset by properties relevant to your target,
   such as molecular size, lipophilicity, charge, reactivity and availability.
   ZINC's [tranche guide](https://wiki.docking.org/index.php/Tranche_Browser)
   explains these choices. Start with a small selection to see its file count
   and storage needs before expanding the screen.
2. Choose **AutoDock PDBQT** (`.pdbqt.gz`) as the download format and **WGET**
   as the method. Download the generated list of commands. If that export is
   already one `wget` command per line for `.pdbqt.gz` archives, copy or rename
   it to `docking/ligands.wget`, replacing the repository's one-archive example.
   Run this from an editable checkout root; for a regular installation, use
   the workspace path printed by `ultidock doctor` before `docking/`:

   ```bash
   cp /path/to/zinc-download-list docking/ligands.wget
   ```

3. Place your receptor as in step 2, then run `ultidock run`. Setup executes
   the list, stores archives in `LIGANDS_DIR`, and the pipeline extracts and
   docks the resulting PDBQT ligands. The normal `auto` hardware mode uses
   available GPUs, including multiple detected NVIDIA GPUs.

The manifest is a list of shell commands, so check the export before using it:
keep only `wget` lines for the intended `.pdbqt.gz` files, with the syntax
shown above. A URL-only list or
an export for SMILES, SDF, MOL2 or DOCK DB2 is not this direct input route.
The [ZINC20 3D archive](https://cache.docking.org/3D/) also exposes available
PDBQT tranche files if you need to inspect a source URL.

### Use prepared local ligands

For local ligands in an editable checkout:

```bash
mkdir -p docking/LIGANDS_DIR
cp /path/to/ligands/*.pdbqt docking/LIGANDS_DIR/
```

The current setup runs the bundled download list even when local ligands
exist. For a local-only run, pass `--skip-wget`; it skips the download, not
docking. Use `--wget FILE` to select a different manifest. Every top-level
PDBQT file directly in `LIGANDS_DIR` is selected, so keep only the intended
ligands there. For ligand chemistry and file-format details, see
[Ligand inputs](preparing-ligands.md).

For a selection from NCBI's [PubChem Advanced Search](https://pubchem.ncbi.nlm.nih.gov/docs/advanced-search),
download 3D structures where available, then prepare ligand PDBQT files as
described in [Ligand inputs](preparing-ligands.md). PubChem's 3D export is SDF,
so renaming that download to `ligands.wget` will not work.

## 4. Run the pipeline

With a receptor and `ligands.wget` in the active workspace:

```bash
ultidock run
```

For local ligand files without the bundled download:

```bash
ultidock run --skip-wget
```

To use a different download manifest explicitly:

```bash
ultidock run --wget /path/to/my-ligands.wget
```

The default `--mode auto` selects a visible NVIDIA CUDA or OpenCL GPU and
falls back to CPU Vina when none is detected. With multiple NVIDIA GPUs,
AutoDock-GPU work is distributed across the detected devices. Use
`--mode gpu` when GPU execution is required; it stops if no compatible GPU is
detected. `--grid-mode centers` is the default CaV-EMPS site search: the
built-in finder proposes docking regions from the receptor.
`--grid-mode blind` uses a broad whole-receptor box instead of a cavity finder.
For fpocket or P2Rank, see [Binding-site discovery](binding-site-discovery/index.md).

The command prints progress through preparation, grids, docking and analysis.
Raw poses (saved 3D ligand coordinates) go to `docking/DOCKING_DIR/`, grids
beside the receptor, and a SQLite score database and CSV tables to
`docking/RESULTS_DIR/` in the active workspace. Use
[Results and reports](results-reports.md) to inspect them.

## Optional: continue into molecular dynamics

On this branch, add `--md-config ./protocol.json` to continue the current screen
through GROMACS system preparation, minimization, NVT and NPT. First install
the [MD tools](../getting-started/md-installation.md) and complete the reviewed
chemistry/system inputs in [Docking to MD](molecular-dynamics/docking-to-md.md).
The protocol names one receptor and a matching docking engine. Production
starts by resuming the printed job after an equilibration review.
