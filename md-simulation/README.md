# Docking to GROMACS

This first CLI integration selects up to five distinct ligand IDs for **one
receptor and one docking engine**. It uses each ligand's lowest-scoring pose
across sites and models, builds a separate all-atom complex, and executes
minimization → restrained NVT → restrained NPT → unrestrained production.
All workflow code lives here. Default runs are in `md-simulation/workspace/`;
regular package installs put that directory in Ultidock's managed home. No GUI
or automatic submission of long production jobs is involved.

## Dependencies

GROMACS, AmberTools (`antechamber`, `parmchk2`, `tleap`, `sqm`), ACPYPE, Open
Babel, RDKit, NumPy and SciPy must be available in the active environment.
Use GROMACS 2021 or newer for C-rescale pressure coupling. Follow
[Molecular Dynamics Installation](../README.md#molecular-dynamics-installation)
in the main README for native GROMACS/Open Babel packages, pip installation
in your existing Python environment, and standalone AmberTools. Conda is an
optional alternative. Check the complete toolchain with `ultidock doctor`.

`python -m pip install -e '.[md]'` installs the Python additions;
the external GROMACS, AmberTools and Open Babel tools are still required.
`md run` accepts `--gmx`, `--acpype` and `--obabel` executable overrides. Use
the thread-MPI `gmx` executable; this runner does not launch MPI jobs.
Dependencies are never silently downloaded by the workflow.
`ultidock doctor` reports docking, receptor preparation and MD dependencies
together, including the active Python environment and both workspaces. Its
`--gmx`, `--acpype` and `--obabel` options can check the same executable
overrides used by `md run`. Missing MD dependencies are reported without
preventing diagnostics for docking-only environments.

## Continue directly from docking

Supply a completed protocol when starting the docking pipeline:

```bash
ultidock run --mode cpu --skip-wget --md-config /path/to/protocol.json
```

The runner checks MD tools and protocol inputs before setup, then checks the
docking engine. Supply chemical inputs for every candidate ligand so the
eventual top-ranked successes can be parameterized. After docking, the handoff
validates the actual receptor, selected poses, ligand chemistry and atom maps,
prepares an MD job and runs through NPT automatically. `--skip-analysis`
does not skip MD.

Screening continues after failed ligand/receptor/site cases. Only successful
outputs enter MD selection; a successful pose at another site remains eligible.
The docking manifest and its `.failures.csv` sibling report skipped cases.
If all docking cases fail, analysis and MD are skipped and the command returns
failure. If the requested MD receptor has no successes, other receptors still
finish and the MD handoff is marked `skipped`. An MD preflight failure also
disables MD without cancelling docking. A failure during MD execution returns
a failure exit status and preserves the job for diagnosis.

All public docking modes (`known-site`, `cavity`, `blind`, `fpocket`, `p2rank`
and `run p2rank/fpocket`) accept `--md-config` (alias `--md`), `--md-through`,
`--md-work-dir`, and `--md-gmx/--md-acpype/--md-obabel`. The default MD endpoint
is `npt`; choose `prepare`, `build`, `em` or `nvt` to stop earlier. Production
requires reviewing and resuming the resulting job.

A docking invocation writes `docking-run-*.json` listing its successful
output artifacts, hashes and failures. Integrated MD ranks only that list, even if
the docking directory contains older results. The sibling `*.md.json` records
the MD job path and handoff status, including a failure reason if the final
preflight rejects the staged receptor before a job is created. The prepared job snapshots the manifest
and binds to the receptor actually staged for docking. MD outputs remain
under this directory's `workspace/`.

For existing results, the same handoff is available directly:

```bash
ultidock md run --config /path/to/protocol.json --docking-dir /path/to/run/docking
```

Add `--pose-manifest /path/to/docking-run-TIMESTAMP.json` to select one recorded
run; otherwise all matching results in that directory are considered. Add
`--work-dir` for a new destination under `md-simulation/`. The CLI prints the
exact job path and resume command, including executable overrides. An MD
execution failure preserves the prepared job for diagnosis/resume; interrupted
system builds still require a fresh job as described below.

## Separate selection and preparation

These commands remain available for inspecting the handoff step by step:

```bash
ultidock md select --docking-dir /path/to/run/docking \
  --receptor-id protein --engine vina
ultidock md prepare --docking-dir /path/to/run/docking \
  --config /path/to/protocol.json
ultidock md run /path/printed/by/prepare --through npt
# After inspecting the built system and equilibration:
ultidock md run /path/printed/by/prepare --through production --equilibration-reviewed
```

Copy and fill in [soluble.json](examples/soluble.json) or
[membrane.json](examples/membrane.json). These are **input templates**, not
validated protocols for the example receptor. Paths are relative to the JSON
file. Fill in the actual lipid counts, atom counts, orientation, hydration and
nonbonded settings; the example numbers do not certify a parameter set.
`--work-dir` can choose a new directory under `md-simulation/`.

The prepared `job.json` preserves scores, binding-site IDs, model/run IDs,
artifact paths and hashes. Source data, selected coordinates and MDP files
are copied into the job. Inputs and workflow code are checked for changes
before execution; changed inputs require a fresh job. Logs and command
arguments are recorded per ligand. A failed MD system stops its own stages;
the runner continues other prepared systems, writes `run_summary.json` and
returns failure after processing the batch if any system failed. Shared input
integrity or dependency failures stop execution before simulation. The runner
never passes `-maxwarn` to GROMACS.

New docking runs record `.md.json` provenance beside each output. AutoDock-GPU
provenance hashes both the score XML and its coordinate DLG, and MD rejects
changed or missing coordinates. Older sidecars without a DLG hash cannot
provide that coordinate-integrity check. The database
`docking_file` field now points to the **output** PDBQT or DLG container; its
`Model` identifier selects the model/run. Older databases used input-ligand
paths, so MD selection does not trust that column. Historical names containing
only `S1__ligand_hash` require `--legacy-single-receptor`: first verify the
directory contains results for exactly the declared receptor. Different
engines are ranked independently. Six pockets for one ligand produce one MD
system, not six candidates. Fewer than five ligands produce fewer systems.

Vina coordinates come from the scored `MODEL`. AutoDock-GPU coordinates come
from the exact XML run ID in its paired DLG and the XML/DLG energies must agree.
A GPU `-best.pdbqt` alone has insufficient run provenance. Missing/truncated
coordinates, non-finite scores and mixed receptor preparations stop selection.

## Ligand chemistry and supported force fields

**Coordinates come from the docking output.** A single-molecule SDF supplies
bond orders, formal charge and explicit stereochemistry, not the ligand's
placement. For each selected ligand provide `net_charge` and `atom_map`:

```json
"atom_map": {"1": 0, "2": 5, "3": 2}
```

Keys are PDBQT **heavy-atom serials**, not line numbers. Values are zero-based
atom indices in the SDF, including any explicit hydrogen indices in its atom
table. The map must cover every heavy atom exactly once. Establish the map
from the original chemical preparation; coordinate proximity and repeated
PDBQT atom names are insufficient. Verify the same protomer/tautomer and
stereoisomer was docked. The workflow checks elements, formal charge,
stereochemistry and plausible mapped bond lengths, then restores missing
hydrogens without moving heavy atoms. Explicitly docked hydrogens are retained
when they agree with the chemical graph. Docking partial charges are not reused as MD charges.

Open Babel converts the posed molecule to MOL2; ACPYPE runs AmberTools to
generate GAFF2/AM1-BCC ligand parameters. Unique atom names allow us to restore
the posed coordinates in ACPYPE's output atom order even if charge calculation
or topology generation reordered or moved atoms. Missing/renamed atoms, wrong
charges or incompatible `[ defaults ]` stop the build. Inspect the AmberTools
logs and parameters; successful typing does not establish accuracy for every
chemical group. This route supports closed-shell organic, noncovalent ligands.
Metals, covalent ligands, radicals, required cofactors and modified residues
need a separate parameterization/system-building provider.

| Protein force field | Allowed water models | Parameters |
| --- | --- | --- |
| `amber99sb-ildn` | `tip3p`, `tip4pew`, `spce` | GROMACS built-in |
| `amber14sb` | `tip3p`, `tip4pew` | Supply `force_field_dir` ending in `amber14sb.ff` |
| `amber19sb` | `opc`, `opc3` | Supply `amber19sb.ff` and matching `water_gro` |

These allowed combinations delimit this provider; they are not a ranking of
water models. An external force-field directory must include the requested
water and ion definitions and is copied into the system. Generic `TIP4P` is
not silently substituted for TIP4P-Ew. CHARMM36m/CGenFF, OPLS and GROMOS are
explicitly rejected by this automatic GAFF2 provider, rather than generating
incompatible topologies. They require future provider implementations. See
[ACPYPE's parameterization pipeline](https://acpype.readthedocs.io/en/latest/)
and [GROMACS force-field guidance](https://manual.gromacs.org/current/user-guide/force-fields.html).

## Receptor and membrane preparation

Provide the exact docking receptor PDBQT and a reviewed, single-model,
protein-only PDB in that same coordinate frame. Heavy-atom identities and
coordinates must agree. Resolve alternate conformations and missing atoms
before docking; changing the binding site afterward requires redocking.
Review protonation, termini and disulfides before setting
`receptor_reviewed=true`. `pdb2gmx_answers` can provide a newline-separated
sequence of answers for a reviewed interactive choice; all choices are logged.
Missing hydrogens are added by `pdb2gmx`; its output topology and chemical
states must be inspected. Required cofactors/structural waters are not silently
discarded by this protein-only builder. See
[pdb2gmx](https://manual.gromacs.org/current/onlinehelp/gmx-pdb2gmx.html).

For membrane proteins such as SERT, use `system_type=membrane`. This provider
requires a **reviewed all-atom, lipid-only periodic bilayer seed** and its
compatible, flattened GROMACS lipid include files; it does not generate or
parameterize arbitrary lipids. The supplied seed may contain multiple species:
describe their contiguous coordinate blocks with molecule names, counts and
atoms per molecule. Atom names/order must match the topology. Label its
protein-force-field compatibility explicitly; GAFF2 ligand compatibility does
not make unrelated lipid parameters compatible.

Supply a rigid 4×4 `transform`, acting on **Angstrom** coordinates, that places
the docking protein into the bilayer's frame. The same rotation/translation is
applied to every docked ligand. The bilayer normal must be Z. Set
`orientation_reviewed=true` after checking insertion depth, orientation,
periodic clearance and the binding pocket. Do not simply center SERT on a
bilayer by its geometric centroid. Membrane seeds currently use orthorhombic
cells, including the rectangular case of `box_shape=triclinic`; tilted cells
are rejected. Dodecahedra are supported for soluble complexes only.

The builder removes **whole** overlapping lipids without shifting the complex,
reports upper/lower leaflet counts, solvates and excludes water in the declared
`hydrophobic_z_nm` slab. `water_pores` supplies optional aqueous pore cylinders
as `{"xy_nm": [x, y], "radius_nm": r}`. These are in **nm** and in the seed
frame. Review pore hydration and conserved waters explicitly for transporters.
Inspect lipid packing, leaflet balance and area per lipid: removing lipids
creates a starting configuration that needs equilibration. The user supplies
`mdp_nonbonded` settings from the lipid parameter authors; the workflow does
not apply a universal dispersion treatment to every lipid force field.

## Salt, equilibration and restart

`salt_molar` specifies added NaCl plus separately added neutralizing ions.
Soluble `salt_basis=box-volume` uses GROMACS `genion -conc`. Membrane default
`water-count` estimates salt pairs from water count using 55.5 M water, avoiding
counting lipid volume as aqueous solvent. Finite box rounding, solvent density
and neutralization affect the realized ionic strength; this is not an exact
aqueous activity specification. The chosen basis, pair count and neutralizing
charge are recorded in `salt.json`. See
[genion's concentration definition](https://manual.gromacs.org/current/onlinehelp/gmx-genion.html).

The baseline uses a timestep ≤2 fs with hydrogen-bond constraints, a single
v-rescale thermostat bath, protein/ligand heavy-atom restraints during NVT/NPT,
and C-rescale pressure coupling (isotropic soluble; semi-isotropic membrane).
NVT generates velocities once; NPT and production carry the preceding
checkpoint. Production removes restraints. Durations are configurable and
**do not establish equilibration by themselves**. Review temperature, density,
pressure trends, ligand pose, contacts, membrane packing, pore hydration and
replicate behavior. This integration does not yet calculate binding free
energies or conclude that a stable trajectory validates a docking score.
See [GROMACS ensemble controls](https://manual.gromacs.org/current/user-guide/mdp-options.html).

`run` defaults to NPT. Production requires `--through production
--equilibration-reviewed`. Each stage writes a status and stops on command
failure, missing outputs or unconverged minimization. Dynamics can resume
using its original TPR/checkpoint; completed outputs are checked before reuse.
A failed/interrupted system build requires a fresh job because build commands
modify topology and solvent in place. Energy traces (`*_energy.xvg`), processed
MDP files, GROMACS logs and `commands.jsonl` are retained for inspection.

## Verification

Run `python -m pytest md-simulation/tests`. Tests cover exact scored pose
extraction, distinct-ligand ranking, original-vs-docked coordinates, chemistry
mapping, incompatible states, topology ordering, salt defaults, whole-lipid
removal, core/pore hydration, ensembles and stop/restart behavior. A small
real GROMACS/AmberTools system and system-specific membrane equilibration are
required before relying on production results; unit tests do not validate the
physical model.

The opt-in `ULTIDOCK_MD_NATIVE_TEST=1 python -m pytest md-simulation/tests/test_workflow.py`
test runs actual ACPYPE/AmberTools and GROMACS through 2 ps each of NVT, NPT and
production on an artificial ethanol/villin complex. Its receptor fixture is
the heavy-atom subset of [RCSB PDB 1VII](https://www.rcsb.org/structure/1VII).
This exercises software and topology consistency, not experimental binding
or membrane physics. It was checked with GROMACS 2025.4 and ACPYPE 2026.9.4.
