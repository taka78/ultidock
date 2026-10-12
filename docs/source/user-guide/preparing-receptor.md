# Automatic receptor preparation

For a normal `ultidock run`, place a receptor `.pdb`, `.mol2` or `.pdbqt`
in the active workspace's `docking/MACRO_MOL_DIR/`. The pipeline prepares and
checks it automatically before building grids. Start with
[Dock your own molecules](start-docking.md); the details below explain what the
preparation stage does and how to inspect a problem input.

For the complete algorithm, where it runs in the pipeline, formatting and
altloc policies, reproducibility conditions and sources, see
[How MolGuard canonicalizes molecular inputs](../scientific-background/molguard-canonicalization.md).

Choose a receptor structure appropriate for your question: chain/assembly,
conformational state, retained cofactors or ions, waters, and protonation matter.
Keep the source and record your edits. Split multiple structural models into
separate receptor files rather than combining coordinates.

For an optional manual check, MolGuard accepts PDB, MOL2 and PDBQT, including
gzip-wrapped inputs:

```bash
molguard receptor prepare source.pdb -o receptor.pdbqt
molguard pdbqt check receptor.pdbqt
```

PDB conversion uses an installed Meeko receptor command when available, with an
Open Babel fallback. MOL2 requires an appropriate converter. A custom command can
be supplied through `--receptor-prepare-command` on the pipeline with `{input}`,
`{output}`, and optional `{seed}` placeholders. Sanitization may resolve altlocs,
remove solvent, rename residues or invoke documented rescue actions. Review the
terminal warnings; do not assume the model is chemically unchanged.

Existing PDBQT files are checked and canonicalized without adding atoms or
recalculating charges. This branch does not recover missing donor hydrogens from
an existing PDBQT. AutoGrid requires HD/HS donor hydrogens: if they are absent,
prepare a hydrogen-complete receptor from its appropriate source PDB. Keep a
copy of the original before in-place preparation; no `.prep.json` recovery
report or automatic recovery backup is produced here.

Canonicalization sorts and renumbers atoms, rewrites representable numeric
fields in AutoDock's fixed-width layout and reports a SHA-256 digest. Its lint
gate rejects numeric-format errors, nonfinite values and missing atom-type
fields. It does not validate model boundaries or torsion trees: only atom
records and a generated remark are written. Supply a single rigid receptor,
and keep ligand/flexible-receptor files out of this path. Unsupported receptor
atom types are checked separately during grid preparation.
Already prepared receptors are not automatically re-protonated. Review
bond-perception and any Meeko residue-deletion rescue warning as a model
change, not a silent format correction.

For pipeline runs, place receptors in `MACRO_MOL_DIR` or pass an absolute
`--macro-mol-dir`. Each filename stem identifies its receptor and site folders.
An existing same-stem PDBQT is preferred over raw inputs unless
`--force-receptor-prep` is supplied. `--receptor-prep-mode off` disables preparation,
not the need for valid engine inputs.

The filename stem is the receptor identity in site folders, grid caches,
database records and output filenames; a fixed name such as `receptor.pdb`
is unnecessary. Multiple stems may be processed together. When both raw
PDB and MOL2 exist for one stem, PDB has priority. Gzip inputs retain their
format suffix during conversion. Pipeline preparation records failed inputs;
other valid receptors can continue through grids and docking. The run manifest
and failure CSV identify excluded receptors and cases.

To canonicalize an existing receptor explicitly:

```bash
molguard receptor canonicalize receptor.pdbqt -o receptor_canon.pdbqt
```

For matching output bytes across repeated calls, add `--timestamp FIXED`.
The direct command otherwise inserts the current UTC time into its remark,
which changes the output hash. This hash identifies file bytes, not chemical
equivalence or biological suitability.

Preserve the raw structure, prepared PDBQT and conversion warnings with the
result. For MD, additionally provide a reviewed protein PDB whose heavy atoms
match the prepared docking receptor in the same frame. Read
[Soluble MD systems](molecular-dynamics/soluble-systems.md) before constructing
that protocol; automatic docking preparation can make structural changes that
must be reconciled with the MD input.

Keep a single reviewed receptor model with valid coordinates, charges and
engine-compatible atom types. See [Troubleshooting](../reference/troubleshooting.md)
for preparation and grid errors.
