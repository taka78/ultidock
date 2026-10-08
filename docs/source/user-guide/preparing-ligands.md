---
orphan: true
---

# Ligand inputs

The docking runner automatically discovers prepared ligand PDBQT files and
extracts configured archives. It does not turn an
arbitrary SMILES or SDF into a scientifically curated screening library for you.
Choose stereochemistry, protonation/tautomer states, conformers and charge models
before screening, and retain the original molecular identifiers and structures.

## Choose a library source

For a large, prebuilt docking library, use the [ZINC20 3D tranches](https://cache.docking.org/3D/)
or [ZINC-22 tranche browser](https://cartblanche22.docking.org/). Select the
chemical-property range that fits the target, then export AutoDock
`.pdbqt.gz` with `wget` commands. Save a list containing only those commands
as `docking/ligands.wget`; the [docking guide](start-docking.md)
shows the short workflow. Ultidock downloads the selected archives and splits
their prepared PDBQT ligands.

NCBI's [PubChem Advanced Search](https://pubchem.ncbi.nlm.nih.gov/docs/advanced-search)
can also define a compound subset. Its [3D structure downloads](https://pubchem.ncbi.nlm.nih.gov/docs/3d-structure-viewer)
are SDF rather than PDBQT; prepare and convert them before putting the resulting
PDBQT files in `LIGANDS_DIR`. A PubChem SDF export cannot be used as
`ligands.wget` by renaming it.

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
that directory; an unrelated molecule left there will also be docked. When
local ligands exist, setup still runs the bundled `ligands.wget` download list.
Pass `--skip-wget` to use only the local files. Compressed archives and `vina_split`
are supported by the extraction stage, but verify the resulting filenames and
ligand count before a large run.

For archive downloads, put one `wget` command per line in
`docking/ligands.wget` or pass `--wget FILE`. A bare URL is not a valid line:
setup executes each nonempty line as a shell command and directs its output to
`LIGANDS_DIR`. If an exported command uses `-O`, put the URL before `-O`
and the output filename, as in `wget URL -O archive.pdbqt.gz`; the current
downloader parses that order. Keep the downloads flat in `LIGANDS_DIR`, since
the extraction stage scans its top level. The repository includes a
one-archive example in the default
manifest. If a download fails, check the source URL and network access or
stage a verified archive locally. The extraction stage expands archives, calls
`vina_split` and removes the temporary unsplit `.pdbqt` so it cannot be
docked as an extra ligand. The original `.gz` archive is retained by default
for ordinary runs. Benchmark mode uses lean artifact retention;
`--no-keep-artifacts` selects that behavior for an ordinary run.
Check the number and names of split ligands against your manifest.

The [D2 example](../tutorials/small-molecule-docking.md) explicitly stages three
ligands from its manifest. Haloperidol's supplied source is PubChem CID 3559;
escitalopram and morphine are comparison molecules, not validated decoy labels.
