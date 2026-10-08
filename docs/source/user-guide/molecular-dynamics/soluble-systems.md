# Soluble systems and ligand chemistry

Use `system_type="soluble"` for a reviewed protein–ligand complex simulated in
water. Start from this {download}`soluble protocol template <../../../../md-simulation/examples/soluble.json>`.
It contains placeholders and unset review flags; fill them for your system
before following [Docking to MD](docking-to-md.md).

```{literalinclude} ../../../../md-simulation/examples/soluble.json
:language: json
```

All file paths resolve relative to the protocol JSON. The full field meanings
and defaults are in the [protocol reference](../../reference/md-protocol.md).

## Match the receptor used for docking

`receptor_id` is the receptor PDBQT filename stem. Supply that exact PDBQT in
`docking_receptor_pdbqt`, plus a reviewed, single-model **protein-only** PDB in
`receptor_pdb`. Heavy-atom identities and coordinates must agree in the same
coordinate frame. The integrated handoff also checks the receptor actually
staged for the docking invocation.

Resolve alternate conformations and missing heavy atoms before docking. Review
protonation, termini and disulfides before setting `receptor_reviewed=true`.
GROMACS `pdb2gmx` adds protein hydrogens and creates the topology; inspect its
choices and output. Optional `pdb2gmx_answers` supplies a newline-separated
answer sequence for reviewed interactive choices, and is logged with the command.
See the [GROMACS pdb2gmx guide](https://manual.gromacs.org/current/onlinehelp/gmx-pdb2gmx.html).

Required cofactors, structural waters, metals and modified residues need an
appropriate system-building provider. This protein-only builder does not silently
discard them to obtain a simulation. Changing the binding-site heavy atoms or
their frame after docking requires a matching new docking preparation and run.

## Supply chemistry for the selected ligands

The `ligands` object is keyed by the **docking ligand filename stem**. Each entry
contains a single-molecule `sdf`, an integer `net_charge` and an explicit
`atom_map`. The SDF supplies bond orders, formal charge and stereochemistry.
The coordinates used for MD come from the selected docking pose.

For example, a map fragment looks like this:

```json
"atom_map": {"1": 0, "2": 5, "3": 2}
```

Keys are PDBQT **heavy-atom serials**, not line numbers. Values are zero-based
atom indices in the SDF atom table, including explicit hydrogen entries in that
table. Provide one mapping for every heavy atom, with no duplicates. Derive it
from the original ligand preparation; repeated atom names or coordinate proximity
do not establish chemical identity. The three-entry fragment above is not a
complete map unless the molecule has exactly those three heavy atoms.

Use the same protomer, tautomer and stereoisomer that was docked. Ultidock checks
elements, formal charge, stereochemistry and mapped bond lengths, then restores
missing hydrogens without moving heavy atoms. Explicitly docked hydrogens are
retained when consistent with the graph.

## Generate and inspect ligand parameters

Open Babel converts the chemically complete posed molecule to MOL2. ACPYPE
runs AmberTools for **GAFF2 atom types and AM1-BCC charges**, then creates
GROMACS ligand parameters. Docking partial charges are not reused as MD charges.
Ultidock restores the original posed coordinates in the topology's atom order
after parameterization. Missing or renamed atoms, an incorrect charge, or
incompatible topology defaults stop that system's build.

Inspect the ACPYPE/AmberTools logs and the resulting `ligand.itp` and
`ligand_atomtypes.itp`. The provider supports closed-shell organic, noncovalent
ligands. Metals, covalent attachment, radicals and required cofactors need
another compatible parameterization route.

## Supported protein force fields and water models

| `force_field` | Allowed `water_model` | Extra input |
| --- | --- | --- |
| `amber99sb-ildn` | `tip3p`, `tip4pew`, `spce` | Built-in GROMACS parameters |
| `amber14sb` | `tip3p`, `tip4pew` | `force_field_dir` ending in `amber14sb.ff` |
| `amber19sb` | `opc`, `opc3` | `force_field_dir` ending in `amber19sb.ff`, plus matching `water_gro` |

An external force-field directory must include the requested water and ion
definitions. It is snapshotted into the job; parameter directories containing
symlinks are rejected. The force-field directory basename must match the chosen
name. TIP4P-Ew requires its own definitions; generic TIP4P is not substituted.
CHARMM/CGenFF, OPLS and GROMOS combinations are rejected by this automatic GAFF2
provider. These allowed combinations describe its support, not a ranking of
force fields or water models.

## Solvate, add ions and equilibrate

The builder combines protein and ligand, creates the periodic box, solvates it
and adds NaCl plus neutralizing ions. `padding_nm` must be at least 1.2 nm for
the default nonbonded cutoff. The default soluble `salt_basis="box-volume"`
uses GROMACS `genion -conc`; `salt.json` records the requested concentration,
basis and neutralizing charge. See the [MD methodology](../../scientific-background/molecular-dynamics.md)
for salt and ensemble assumptions.

The template's durations and temperature are starting settings to edit for
your system. Follow [Run stages, results and restart](running-and-results.md)
to check the built complex and equilibration before production.
