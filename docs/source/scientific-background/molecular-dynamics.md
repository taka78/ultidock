# Molecular dynamics after docking

Docking proposes poses and ranks them with an engine-specific score. Molecular
dynamics (MD) follows the motion of a prepared protein–ligand system under a
force field. The `gmx-dev` workflow connects these stages by taking coordinates
from an identified, scored docking pose and preserving their relationship to
the receptor. It does not infer a chemically complete simulation system from
a PDBQT file alone.

## Chemical identity and coordinate integrity

The reviewed protein PDB supplies the protein for GROMACS. Its heavy atoms must
match the docking receptor in the same coordinate frame. Each ligand SDF supplies
bond orders, stereochemistry and explicit chemistry; its atom map links docking
PDBQT serial numbers to SDF atom indices. The selected pose supplies the ligand's
coordinates. Parameterization uses GAFF2 and AM1-BCC through ACPYPE/AmberTools;
the prepared coordinates are restored to the selected pose afterward.

This division lets the workflow check identity and geometry without treating
docking atom types or partial charges as a complete MD topology. Protonation,
tautomer choice, missing residues, cofactors and biological assembly still need
review before the protocol is completed. See [Soluble systems](../user-guide/molecular-dynamics/soluble-systems.md).

## What the stages do

| Stage | Purpose in this workflow |
| --- | --- |
| Build | Generate the protein and ligand topology, combine coordinates, solvate and add ions |
| Energy minimization | Reduce unfavorable contacts before dynamics |
| NVT | Establish the requested temperature at fixed volume, with positional restraints |
| NPT | Relax pressure and volume with positional restraints |
| Production | Continue unrestrained dynamics after an explicit equilibration review |

The generated dynamics use a timestep no greater than 2 fs and constrain
covalent bonds involving hydrogen. NVT initializes velocities; later stages continue from preceding
state/checkpoint files. Temperature coupling uses V-rescale. NPT and production
use C-rescale pressure coupling, isotropic for soluble systems and semi-isotropic
for membrane systems. GROMACS documents the meaning and settings of these
[thermostat and barostat options](https://manual.gromacs.org/current/user-guide/mdp-options.html).

Durations in the templates are starting settings. Completion of NPT is a software
milestone; examine temperature, density, pressure behavior, restraints, contacts
and the geometry of your system before requesting production. The runner records
energy samples and logs, but does not determine scientific equilibration for you.

## Salt and membrane assumptions

Soluble systems request NaCl concentration from GROMACS using the simulation-cell
volume, with charge neutralization added separately. GROMACS describes this
behavior in its [genion reference](https://manual.gromacs.org/current/onlinehelp/gmx-genion.html).
Membrane systems instead estimate salt pairs from the retained water count using
the protocol concentration and a 55.5 M water reference. Neutralizing ions are
additional in both cases; the membrane builder records the resulting counts in
`salt.json`.

A membrane needs compatible, reviewed lipid parameters and a bilayer seed.
The builder removes overlapping whole lipid molecules and excludes core water
outside declared pores. These operations can alter lipid counts and local
packing. Inspect both leaflets, composition, pore hydration and protein
orientation, then equilibrate for the system being studied. The
[membrane guide](../user-guide/molecular-dynamics/membrane-systems.md) explains
the required coordinates and review fields.

## Interpreting the outcome

Pose selection uses the lowest binding-energy pose per distinct ligand, for the
protocol's single receptor and scoring engine. It reads raw scored outputs;
analysis CSV filters do not determine MD selection. This is a reproducible
handoff rule, rather than evidence that the selected pose is biologically correct.

An intact complex during a short trajectory does not establish affinity or
binding stability. This branch does not automatically calculate binding free
energy, compare replicate uncertainty or produce a validated MD ranking. Retain
the protocol, selected pose identities, input hashes and stage outputs when
analyzing trajectories. The [recorded acceptance tests](../development/md-validation.md)
demonstrate software execution and restart integrity with short simulations.
