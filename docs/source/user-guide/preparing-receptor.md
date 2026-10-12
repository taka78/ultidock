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

Existing PDBQT files are normalized, checked and canonicalized. If no HD/HS donor
hydrogens exist, the recovery path uses Open Babel to add hydrogens and recalculate
atom types/charges. A candidate must preserve heavy-atom identities and coordinates.
The original is backed up under `.molguard-backups/`; `<stem>.prep.json` records
recovery and converter diagnostics. This does not repair missing heavy atoms or
prove that a partially hydrogenated structure is complete.

Canonicalization sorts and renumbers atoms, rewrites representable numeric
fields in AutoDock's fixed-width layout and reports a SHA-256 digest. Its direct
command lints formatting first, then retains atom records and a generated remark.
The shared `receptor prepare` path adds numeric repair and rigid-receptor
validation before canonicalization can discard context. That preparation path
rejects missing coordinates/charges, nonfinite values, unsupported receptor
types, multiple models and ligand/flexible-receptor torsion records.
Donor-hydrogen recovery infers bonds from the existing heavy-atom coordinates,
adds hydrogens and recalculates Gasteiger charges. It must preserve each heavy
atom's identity, element and coordinates before replacing the receptor.
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
format suffix during conversion. Preparation collects failures and stops
before grid generation, avoiding a partially prepared receptor queue.

To canonicalize an existing receptor explicitly:

```bash
molguard receptor canonicalize receptor.pdbqt -o receptor_canon.pdbqt
```

For matching output bytes across repeated calls, add `--timestamp FIXED`.
The direct command otherwise inserts the current UTC time into its remark,
which changes the output hash. This hash identifies file bytes, not chemical
equivalence or biological suitability. Use `molguard receptor prepare` when
numeric repair, additional rigid-receptor validation or hydrogen recovery is
needed; the direct canonicalizer does not run those extra steps.

Preserve the raw structure, backup, prepared PDBQT, `.prep.json` and any
conversion warnings with the result.

Unknown atom types, non-finite values, missing required coordinates/charges, and
incompatible torsion/model records fail validation. See [Troubleshooting](../reference/troubleshooting.md).
