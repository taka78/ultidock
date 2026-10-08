# Molecular dynamics with GROMACS

On `gmx-dev`, Ultidock can continue a docking screen into all-atom molecular
dynamics (MD). MD follows how a prepared protein–ligand complex moves over
time under a chosen force field, solvent and simulation protocol.

The workflow selects the lowest-scoring pose of each distinct ligand for
**one receptor and one docking engine**. By default it takes five ligands,
or fewer when fewer qualify. A ligand with poses in several pockets still
produces one MD system. The protocol's `top` field changes the requested count.
Each selected complex is built and run separately.

## Choose a starting point

| Your starting point | Guide |
| --- | --- |
| MD programs are not installed yet | [Install the MD tools](../../getting-started/md-installation.md) |
| Start docking and automatically continue into MD | [Docking to MD, step by step](docking-to-md.md) |
| Prepare a protein–ligand complex in water | [Soluble systems and ligand chemistry](soluble-systems.md) |
| Prepare a membrane protein | [Membrane systems](membrane-systems.md) |
| Review, resume or inspect a prepared job | [Run stages, results and restart](running-and-results.md) |
| Look up JSON fields or command options | [Protocol reference](../../reference/md-protocol.md) and [MD CLI reference](../../reference/md-cli.md) |

## What you supply

Supply a reviewed protein PDB, the exact receptor PDBQT used for docking, and
a single-molecule SDF plus charge and atom map for each ligand that may be
selected. The SDF defines the chemical graph; the actual ligand coordinates
come from the scored docking output. A membrane system additionally needs a
reviewed bilayer, compatible lipid parameters and a common coordinate transform.
The JSON [protocol](../../reference/md-protocol.md) connects these inputs.

Ultidock validates and snapshots the inputs, restores ligand chemistry,
generates GAFF2/AM1-BCC ligand parameters, and builds a GROMACS complex. It
runs energy minimization (EM), restrained equilibration at constant temperature
and volume (NVT), then restrained equilibration at constant temperature and
pressure (NPT). A new job stops at NPT by default. After reviewing equilibration,
resume it to run unrestrained production dynamics.

Docking's worker pool can screen many receptors and ligands in one run.
The current MD handoff targets the one `receptor_id` and `engine` in its
protocol, and executes prepared complexes sequentially. See
[HPC / batch screening](../hpc-batch-screening.md) for the screening workflow.
The supplied example protocols are editable input templates; the bundled
docking examples do not include complete reviewed MD systems.

```{toctree}
:maxdepth: 1

docking-to-md
soluble-systems
membrane-systems
running-and-results
```
