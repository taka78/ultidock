# Preparing a receptor

Choose a receptor structure appropriate for your question: chain/assembly,
conformational state, retained cofactors or ions, waters, and protonation matter.
Keep the source and record your edits. Split multiple structural models into
separate receptor files rather than combining coordinates.

MolGuard accepts PDB, MOL2 and PDBQT, including gzip-wrapped inputs:

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

For pipeline runs, place receptors in `MACRO_MOL_DIR` or pass an absolute
`--macro-mol-dir`. Each filename stem identifies its receptor and site folders.
An existing same-stem PDBQT is preferred over raw inputs unless
`--force-receptor-prep` is supplied. `--receptor-prep-mode off` disables preparation,
not the need for valid engine inputs.

Unknown atom types, non-finite values, missing required coordinates/charges, and
incompatible torsion/model records fail validation. See [Troubleshooting](../reference/troubleshooting.md).
