# Membrane regression and native tests

Run these from a source checkout in an environment with Ultidock's development
and MD dependencies. Ordinary regression tests require no downloaded structures:

```bash
python -m pytest -q md-simulation/tests/test_membrane.py md-simulation/tests/test_systems.py
```

They exercise periodic lipid overlap removal, preservation of whole molecules,
both retained leaflets, periodic pore-water selection, three- and four-site
water grouping, topology/coordinate ordering, box geometry and periodic padding.

## Native fixture

The opt-in test combines the authors' **Slipids 2016 POPC bilayer at 303 K**
with the **OPM-oriented glycophorin A dimer, 1AFO**. A synthetic benzene pose
exercises ligand parameterization and preservation of the common coordinate
frame. It is not an experimentally supported binding pose or a fresh docking
run. The ordinary docking-to-MD handoff is covered separately in
[the acceptance report](ACCEPTANCE.md).

Sources and rationale:

- [Slipids authors' dataset](https://zenodo.org/records/1149623): force-field
  parameters, POPC topology and equilibrated bilayer coordinates. The supplied
  POPC topology explicitly lists Amber99SB-ILDN compatibility. Its defaults
  match Amber/GAFF2: combination rule 2, LJ scaling 0.5 and electrostatic scaling
  0.8333. See also [the force-field paper](https://doi.org/10.1021/ct300342n).
- [OPM](https://opm.phar.umich.edu/) supplies the orientation of the
  [glycophorin A transmembrane dimer](https://www.rcsb.org/structure/1AFO).
  The downloaded structure places the membrane midplane at Z=0 with a
  1.595 nm half-thickness. Dummy membrane markers are excluded from the protein.
- The test uses PME, a 1.4 nm LJ cutoff, potential shift and energy/pressure
  dispersion correction. Both real-space cutoffs are 1.4 nm for Verlet.
  The authors discuss the LJ cutoff and dispersion correction in
  [their parameter guidance](https://www.mail-archive.com/gmx-users@gromacs.org/msg57875.html).

Download only the small source archive and oriented protein, then prepare the
fixture. Downloads happen explicitly here; pytest never downloads data:

```bash
mkdir -p md-simulation/workspace/membrane-sources
curl -fL https://zenodo.org/records/1149623/files/Slipids_2016.tar.gz \
  -o md-simulation/workspace/membrane-sources/Slipids_2016.tar.gz
curl -fL https://opm-assets.storage.googleapis.com/pdb/1afo.pdb \
  -o md-simulation/workspace/membrane-sources/1afo-opm.pdb
python md-simulation/tests/prepare_membrane_fixture.py \
  --slipids-archive md-simulation/workspace/membrane-sources/Slipids_2016.tar.gz \
  --opm-pdb md-simulation/workspace/membrane-sources/1afo-opm.pdb \
  --output md-simulation/workspace/membrane-fixture
ULTIDOCK_MD_NATIVE_TEST=1 \
ULTIDOCK_MD_MEMBRANE_FIXTURE="$PWD/md-simulation/workspace/membrane-fixture" \
  python -m pytest -q -s md-simulation/tests/test_membrane_native.py
```

The helper requires a new output directory and checks pinned source checksums.
An upstream file change requires review, not updating the checksum blindly.
`review.json` records sources, transformations, assumptions and file hashes;
the native test verifies those hashes before using the fixture.

The helper preserves the original XY lattice, removes source water and extends
only the aqueous Z space to 11 nm. It retains the OPM protein orientation and
applies a common translation to protein and ligand. It prefixes lipid atom-type
names to prevent collisions with protein/ligand types, excludes parameter terms
for unrelated lipid types, and preserves all applicable POPC terms and numerical
values. This avoids an unused duplicate torsion in the complete source bundle
without bypassing GROMACS warnings.

The protein contains both chains, residues 66–101. Charged fragment termini
are used; pdb2gmx selected neutral epsilon-protonated HIS66/HIS67 in both chains
in the recorded run, giving protein charge +4. The ligand has charge zero.
The native test checks salt pairs from the water count, additional neutralizing
ions, transformed ligand coordinates, whole retained lipids, absence of core
water outside declared pores, semi-isotropic NPT, finite energy samples and
absence of LINCS warnings. It runs actual GROMACS and ACPYPE/AmberTools through
build, minimization, 2 ps NVT, 2 ps NPT and 2 ps production.

## Scientific scope

Review flags in this fixture apply to this documented **software test**. They
do not approve its geometry, fragment termini, histidine choices or short
equilibration for research production. Lipid deletion leaves a packing defect
that needs substantially more equilibration and inspection. The test does not
validate a physiological glycophorin–benzene interaction, SERT orientation,
transporter pore hydration, mixed-lipid behavior or binding stability. A SERT
production system still needs its own reviewed biological assembly, chemical
states, membrane composition/orientation and equilibration protocol.
