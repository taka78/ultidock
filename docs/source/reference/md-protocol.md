# MD protocol reference

MD takes a JSON protocol through `--config` or the docking option
`--md-config`. Paths in that JSON resolve relative to the JSON's own directory.
After preparation, `job.json` snapshots the resolved inputs and stage settings.
Use [Docking to MD](../user-guide/molecular-dynamics/docking-to-md.md) for the
command sequence.

Download the {download}`soluble template <../../../md-simulation/examples/soluble.json>`
or {download}`membrane template <../../../md-simulation/examples/membrane.json>`.
Their review flags, atom maps and scientific settings require completion.

## Identity and molecular inputs

| Field | Required value or behavior |
| --- | --- |
| `system_type` | Explicit `soluble` or `membrane` |
| `receptor_id` | Receptor filename stem in docking results |
| `engine` | `vina` or `adgpu`; choose one scoring engine |
| `receptor_pdb` | Reviewed single-model, protein-only PDB |
| `docking_receptor_pdbqt` | Exact docking preparation in the same coordinate frame |
| `receptor_reviewed` | Must be `true` after completeness, protonation, termini, disulfides and cofactors have been reviewed |
| `ligands` | Object keyed by ligand filename stem; each entry has `sdf`, `net_charge` and `atom_map` |
| `ligands.<id>.sdf` | Single-molecule SDF defining bond orders, formal charge and stereochemistry |
| `ligands.<id>.net_charge` | Integer molecular charge |
| `ligands.<id>.atom_map` | Complete one-to-one mapping: PDBQT heavy-atom serial string → zero-based SDF atom index |
| `pdb2gmx_answers` | Optional newline-separated answers for reviewed GROMACS protein choices |

Per-ligand chemistry is checked when candidates are checked or selected.
Invalid selected candidates are reported and skipped individually. Shared
receptor, protocol and integrity errors stop the MD preparation. See
[Soluble systems and ligand chemistry](../user-guide/molecular-dynamics/soluble-systems.md).

## Force field and solvent

| Field | Default and restrictions |
| --- | --- |
| `force_field` | `amber99sb-ildn`; also supports `amber14sb` or `amber19sb` with external parameters |
| `water_model` | `tip3p`; allowed combinations are listed below |
| `force_field_dir` | Required for Amber14SB/Amber19SB; directory name must be `<force_field>.ff` and include water/ion definitions |
| `water_gro` | Optional matching solvent coordinate template; required for OPC/OPC3 |
| `box_shape` | `dodecahedron` if omitted; allowed: `cubic`, `dodecahedron`, `triclinic`; membrane must use a planar orthorhombic cell |
| `padding_nm` | 1.2 nm; minimum 1.2 nm and at least the largest membrane cutoff |
| `salt_molar` | 0.15 M added NaCl; finite and nonnegative; neutralizing ions are added separately |
| `salt_basis` | `box-volume` for soluble, `water-count` for membrane; either is accepted |

| Protein force field | Allowed water models |
| --- | --- |
| `amber99sb-ildn` | `tip3p`, `tip4pew`, `spce` |
| `amber14sb` | `tip3p`, `tip4pew` |
| `amber19sb` | `opc`, `opc3` |

The automatic ligand provider uses GAFF2/AM1-BCC. CHARMM/CGenFF, OPLS and
GROMOS are unsupported by this provider. External parameter directories are
copied into the job and must not contain symlinks.

## Sampling and resources

| Field | Default when omitted | Validation or meaning |
| --- | --- | --- |
| `top` | 5 | Positive 32-bit integer; number of distinct ligand IDs requested, not pockets or engine models |
| `threads` | 1 | Positive 32-bit integer; CPU thread budget for each sequential GROMACS simulation |
| `temperature_k` | 310 K | Positive finite target temperature |
| `nvt_ps` | 100 ps | Restrained NVT duration |
| `npt_ps` | 1,000 ps | Restrained NPT duration |
| `production_ns` | 100 ns | Unrestrained production duration after review |
| `dt_ps` | 0.002 ps | Positive timestep, at most 0.002 ps (2 fs) |
| `seed` | 2026 | Positive 32-bit integer used for velocities and ion placement |

Durations must be positive finite values and exact integer multiples of the
timestep after conversion to ps. Template values may differ from omitted-field
defaults: for example, the membrane template supplies longer equilibration and
both templates explicitly request four threads. Durations are user choices,
not an equilibration certificate. [Run stages and results](../user-guide/molecular-dynamics/running-and-results.md)
explains the generated MDPs and resource behavior.

## Membrane object

All membrane input geometry must refer to the supplied seed's frame. See
[Membrane systems](../user-guide/molecular-dynamics/membrane-systems.md) before
setting review flags or using the example values.

| `membrane` field | Meaning |
| --- | --- |
| `gro` | Reviewed all-atom, lipid-only periodic bilayer seed |
| `topology_dir` | Directory containing compatible flattened lipid includes |
| `include_files` | Nonempty ordered list of parameter/molecule include files; nested `#include` directives are rejected |
| `molecules` | Ordered contiguous coordinate blocks: `name`, `count`, `atom_count` per molecule |
| `compatible_force_field` | Must equal the chosen protein force-field name after compatibility review |
| `orientation_reviewed` | Must be `true` after reviewing insertion depth, orientation, periodic clearance and pocket placement |
| `transform` | Finite rigid 4×4 transform on **Å coordinates**, applied to protein and each selected ligand |
| `hydrophobic_z_nm` | Finite lower/upper Z slab bounds in **nm**, in the bilayer frame |
| `water_pores` | Optional cylinders with `xy_nm: [x, y]` and positive `radius_nm`; periodic XY distances are used |
| `clash_distance_nm` | 0.2 nm if omitted; positive lipid/complex overlap threshold |
| `mdp_nonbonded` | Lipid-author settings listed below; required |

The membrane normal is Z. Tilted cells and dodecahedral membrane boxes are
unsupported. Both leaflets must retain lipids. The builder sets water-site
grouping from the selected water model; supply complete seed/topology blocks.

`mdp_nonbonded` requires `rvdw` and `rcoulomb` in nm, `vdw-modifier` and
`DispCorr`. Accepted modifiers are `Potential-shift`, `Force-switch` and
`Potential-switch`; accepted dispersion corrections are `no`, `Ener` and
`EnerPres`. Switching requires `0 < rvdw-switch < rvdw`. Other keys in this
object are rejected. Use values from the actual lipid parameter authors.

## Generated stage parameters

These are generated workflow settings rather than extra JSON fields:

| Setting | Current baseline |
| --- | --- |
| Electrostatics | PME for minimization/dynamics, Verlet cutoff scheme |
| Soluble cutoffs | 1.2 nm Coulomb/Lennard-Jones, potential shift and `EnerPres` dispersion correction |
| Constraints | Bonds involving hydrogen (`h-bonds`), LINCS order 4 |
| Minimization | Steepest descent, up to 50,000 steps, 1,000 kJ mol⁻¹ nm⁻¹ force threshold |
| Temperature | One `System` V-rescale bath, coupling time 0.5 ps |
| Pressure | C-rescale, coupling time 5 ps, 1 bar; isotropic soluble or semi-isotropic membrane |
| Restraints | Protein and ligand heavy atoms during NVT/NPT; removed in production |
| Logging and energy | Every 500 dynamics steps |
| Compressed trajectory | Every 5,000 dynamics steps; at the default timestep this is 10 ps |

NVT generates velocities with the configured seed. Later stages use the preceding
checkpoint. Membrane `mdp_nonbonded` replaces the baseline cutoff/dispersion
settings. Very short smoke runs may be shorter than a trajectory output interval.
Inspect the generated and processed MDP files retained by
[the runner](../user-guide/molecular-dynamics/running-and-results.md).
