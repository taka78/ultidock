# Membrane systems

Use `system_type="membrane"` for a receptor that needs a lipid bilayer. The
workflow combines the docked complex with a **reviewed all-atom, lipid-only
periodic bilayer seed** and supplied compatible lipid parameters. It does not
generate an arbitrary membrane or infer a protein's insertion orientation.

Start from this {download}`membrane protocol template <../../../../md-simulation/examples/membrane.json>`:

```{literalinclude} ../../../../md-simulation/examples/membrane.json
:language: json
```

The SERT names, lipid counts, slab bounds and identity transform are placeholders.
The review flags start false. Complete the receptor and ligand chemistry from
[Soluble systems](soluble-systems.md), then review these additional inputs.

## 1. Supply a compatible bilayer and topology

Set `membrane.gro` to a lipid-only periodic seed. Set `topology_dir` to the
compatible GROMACS lipid parameters and list their ordered, flattened
`include_files`. The include files must contain all dependencies; nested
`#include` directives are rejected during preparation.

Describe the seed's contiguous coordinate blocks in `molecules`. Each block
needs its topology molecule `name`, molecular `count` and `atom_count` per
molecule. Mixed lipid species are supported when the ordered blocks and
parameters match the coordinates. Atom names and ordering must match the
topology. Set `compatible_force_field` to the chosen protein force field only
after confirming its compatibility with the lipid parameters and GAFF2 ligand.

Use an orthorhombic cell with membrane normal along Z. A rectangular cell
specified as `box_shape="triclinic"` is accepted; tilted cells and soluble
dodecahedral boxes are rejected for this membrane provider.

## 2. Orient the docked complex in the bilayer frame

`membrane.transform` is a rigid 4×4 rotation/translation acting on **Ångström
coordinates**. It places the docking protein into the seed's coordinate frame.
Ultidock applies the same transform to every selected ligand, preserving the
relative protein–ligand placement. The rotation must preserve distances and
handedness, and the last row must be `[0, 0, 0, 1]`.

Inspect insertion depth, orientation, periodic clearance and the pocket before
setting `orientation_reviewed=true`. A geometric-centroid placement alone does
not establish a suitable transporter orientation. A compatible bilayer and
orientation prepared for another protein do not validate the SERT template.

## 3. Specify hydration and nonbonded settings

`hydrophobic_z_nm` gives the lower and upper Z bounds of the hydrophobic slab
in **nm**, in the seed frame. The builder excludes complete water molecules in
that slab except inside optional `water_pores`. Each pore is an aqueous cylinder:

```json
{"xy_nm": [3.0, 4.0], "radius_nm": 0.5}
```

Pore locations and radii are in nm and use periodic distances in the membrane
plane. Review pore hydration and conserved waters for transporters. Water
grouping preserves complete three-site or four-site molecules according to
the selected water model.

Set `mdp_nonbonded` using the lipid parameter authors' guidance. Required keys
are `rvdw`, `rcoulomb`, `vdw-modifier` and `DispCorr`; switching modifiers also
need `rvdw-switch`. Cutoffs are in nm. `padding_nm` must cover the largest
cutoff. The example values do not certify a lipid parameter set. See the
[protocol reference](../../reference/md-protocol.md) for accepted values.

## 4. Build and inspect the membrane complex

Follow the [docking-to-MD command sequence](docking-to-md.md) with the completed
membrane protocol. `--md-through build` stops after system assembly and ion
addition when you want to inspect those outputs first.

The builder removes **whole overlapping lipid molecules** using
`clash_distance_nm`, without shifting the complex. It writes retained species
and leaflet counts to `membrane_composition.json`; both leaflets must retain
lipids. It solvates the system, applies slab/pore water selection and adds salt.
The default membrane `salt_basis="water-count"` estimates NaCl pairs from
water count using 55.5 M water, then adds neutralizing ions separately.
The estimate and charge are recorded in `salt.json`.

Inspect lipid packing, leaflet balance, area per lipid, protein/ligand
placement and pore water before continuing. Removing overlapping lipids creates
a starting configuration that requires equilibration. NPT and production use
semi-isotropic pressure coupling so the membrane plane and normal can respond
separately. Review those trajectories before requesting production.

The [MD validation guide](../../development/md-validation.md) documents the
reproducible POPC/glycophorin software fixture. It exercises this builder and
native GROMACS execution; system-specific membrane preparation and equilibration
remain necessary for a research run.
