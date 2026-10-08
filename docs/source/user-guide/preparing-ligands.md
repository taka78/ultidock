# Preparing ligands

The docking runner consumes prepared ligand PDBQT files. It does not turn an
arbitrary SMILES or SDF into a scientifically curated screening library for you.
Choose stereochemistry, protonation/tautomer states, conformers and charge models
before screening, and retain the original molecular identifiers and structures.

Use a ligand preparation tool such as Meeko or AutoDockTools appropriate to your
input and installed version. Consult the tool's help for conversion options.
Treat alternative chemical states as distinct inputs with traceable names.
PDBQT contains atom types, charges and a torsion tree; it is not a lossless archive
of the source molecule's chemistry.

For Meeko in an activated virtual environment:

```bash
python -m pip install meeko
mk_prepare_ligand.py --help
# If that executable is absent, try:
python -m meeko.cli.mk_prepare_ligand -h
```

Prepare source chemistry and conformers before conversion; the exact command
depends on the source format and Meeko version. The main Ultidock run accepts
the resulting PDBQT files, not raw SMILES or SDF.

Check a prepared file before staging it:

```bash
molguard pdbqt check ligand.pdbqt
molguard pdbqt normalize ligand.pdbqt -o ligand-normalized.pdbqt
```

Normalization preserves atom ordering and ROOT/BRANCH/TORSDOF records. It is a
formatting operation, not a way to add missing chemistry. Resolve non-finite
values, missing types, and invalid preparation with the source tool.

Put only the desired `.pdbqt` files in the run's `LIGANDS_DIR`. The pipeline scans
that directory; an unrelated molecule left there will also be docked. Use
`--skip-wget` for a prepared local library. Compressed archives and `vina_split`
are supported by the extraction stage, but verify the resulting filenames and
ligand count before a large run.

The [D2 example](../tutorials/small-molecule-docking.md) explicitly stages three
ligands from its manifest. Haloperidol's supplied source is PubChem CID 3559;
escitalopram and morphine are comparison molecules, not validated decoy labels.
