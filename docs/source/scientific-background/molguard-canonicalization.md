# How MolGuard canonicalizes molecular inputs

MolGuard gives Ultidock a consistent representation of prepared docking inputs.
For a **receptor**, canonicalization selects alternate-location records, orders
atoms, replaces serial numbers and writes a standardized PDBQT. For a **ligand**,
normalization rewrites numeric fields while preserving its atom order and torsion
records. The shared receptor-preparation wrapper adds rigid-receptor validation and,
when needed, donor-hydrogen recovery before final output replacement.

| Operation | Input | Main effect |
| --- | --- | --- |
| `molguard pdbqt check` | Receptor or ligand PDBQT | Report formatting errors and warnings without writing |
| `molguard receptor canonicalize` | Prepared rigid-receptor PDBQT | Select altloc, sort/renumber atoms, reformat and return an output hash |
| `molguard pdbqt normalize` | Prepared ligand PDBQT | Reformat atom records while preserving atom serials/order and other records |
| `molguard receptor prepare` | PDB, MOL2 or PDBQT, including gzip inputs | Normalize/validate prepared receptors, recover absent donor hydrogens when needed, or sanitize/convert raw receptors; finish with canonicalization |

This page describes the implementation on `release/v1.1.2` at commit
`f617684c4de0363f30c54095ca5dc6cd19a99e7d`. The implementation links below
identify the policies used by this version.

## Why these steps are needed

PDBQT extends PDB coordinates with partial charges, AutoDock atom types and
ligand torsion records. Its numeric layout and torsion-tree conventions are
specified in the [AutoDock 4.2.6 User Guide, Appendix I](https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf).
The [wwPDB coordinate specification](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)
defines the underlying coordinate fields, alternate-location identifiers and
insertion codes.

The bundled [AutoGrid parser](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/docking/AUTODOCK_GPU_DIR/autogrid/mainpost1.28.cpp#L875)
reads coordinates, charges and atom types from fixed character offsets. Keeping
those offsets consistent prevents field displacement. MolGuard additionally
checks decimal notation, finite values and ASCII, and standardizes output line
endings. These are local policies. Do not infer that every native parser
rejects every exponent or integer representation in the same way.

The sort order, altloc selection, B-factor handling and timestamp conventions
are Ultidock implementation choices. The external format specifications explain
the fields; they do not prescribe MolGuard's canonicalization algorithm.

## Receptor canonicalization, step by step

The following sequence is implemented by
[`canonicalize_receptor`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/io/pdbqt.py#L462).

### 1. Check the input

`pdbqt_check` decodes the file as ASCII and checks the numeric slices of each
`ATOM`/`HETATM` record. It collects truncated fields, exponent notation, missing
decimal points, unparseable numbers, NaN/Inf and missing atom-type fields.
Short atom records and B-factors above 99.99 produce warnings. A lint error
stops canonicalization or normalization before the output is written.

This checker covers a limited format contract. It does not validate bond orders,
protonation, structural completeness, torsion-tree consistency or membership of
each nonblank atom type in an allowed parameter set. Blank numeric fields can
also fail later when parsed. The shared receptor-preparation path performs
additional validation through `receptor_checks.py` before canonicalization can
remove model/torsion context. Atom-type compatibility is also checked during
[grid preparation](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/docking/make_grids.py#L288).

### 2. Parse atom records and select alternate locations

The canonicalizer parses `ATOM` and `HETATM` records into atom objects. If
nonblank altloc labels exist, it keeps the alphabetically first label **for each
residue**, together with all blank-altloc atoms. Residues are grouped by chain,
sequence number, insertion code and residue name. Other labeled atoms are
discarded with a warning, and retained labels are cleared.

For example, one residue with A/B retains A, while another with only B retains
B. This is a deterministic selection rule, not a highest-occupancy calculation
or a review of biological suitability. Select the intended conformer before
docking when this choice matters. Raw PDB sanitization explicitly prefers A
per residue when available, as described below.

### 3. Sort atoms and replace serial numbers

Atoms are sorted in ascending order by this exact tuple:

```text
(chain_id, residue_sequence_number, residue_name,
 stripped_atom_name, AutoDock_atom_type)
```

The residue number is compared numerically; string fields are compared
lexicographically. This is not a backbone ordering such as N–CA–C–O. Original
serial numbers are excluded, and new serials become `1, 2, …, N`.

The insertion code is retained in the output but is absent from the sort key.
Coordinates are also absent. Equal keys keep their input order through Python's
stable sort, so duplicate keys can prevent arbitrary permutations from producing
identical bytes. This is file canonicalization for a prepared receptor, not a
graph-based molecular identifier or a test of chemical equivalence.

### 4. Rewrite numeric fields

Both receptor and ligand writers use
[`molguard.io.fixedfmt`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/io/fixedfmt.py).
Column numbers below are one-based and inclusive:

| Field | Columns | Width | Decimal places |
| --- | --- | --- | --- |
| X, Y, Z in Å | 31–38, 39–46, 47–54 | 8 each | 3 |
| Occupancy | 55–60 | 6 | 2 |
| B-factor | 61–66 | 6 | 2 |
| Partial charge | 71–76 | 6 | 3 |
| AutoDock atom type | 78–79 | 2 | Text, right-justified |

Formatting uses a decimal point independent of the system locale, excludes
scientific notation and normalizes an exact floating-point `-0.0` to `0.0`.
The writer checks representability instead of allowing numeric fields to grow
into adjacent columns. Rewriting rounds coordinates and charges to the stated
precision; it does not optimize geometry or recalculate charges.

For inputs passing lint, B-factor output is `min(abs(B), 99.99)`. This is a
local metadata policy intended to retain leading space in that field. Preserve
the source structure if its original B-factors are needed for structural analysis.
This cap is not a constraint imposed by the wwPDB specification.

### 5. Write the file and calculate its hash

The output contains a generated `REMARK` with MolGuard's version and stamp,
followed by the sorted atom records. Other source records, including original
remarks, `TER`, `MODEL`, connectivity and torsion records, are not emitted.
Output is ASCII with LF line endings and a final newline. SHA-256 is calculated
over the **complete output bytes** and returned to the caller.

Because this writer collects atom records without validating model boundaries,
provide a single prepared rigid-receptor structure. Separate structural models
before running it. Do not use receptor canonicalization on a ligand or a flexible
receptor file whose topology records need to survive.

## Ligand normalization preserves the torsion tree

[`pdbqt_normalize`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/io/pdbqt.py#L295)
applies the same lint gate and numeric formatter, but writes each record back in
its original order. Atom serials and nonnumeric atom fields are preserved.
`ROOT`, `ENDROOT`, `BRANCH`, `ENDBRANCH`, `TORSDOF` and other non-atom records
are passed through unchanged apart from line-ending normalization.

Ligand branch references depend on atom serials and record placement, which is
why this path avoids sorting or renumbering. A formatting error such as
`EXPONENT` or `NO_DECIMAL` is **rejected** by direct `pdbqt normalize` and
`receptor canonicalize`. Correct the ligand source writer or regenerate a valid
prepared file before normalization. Warning-only B-factor issues can be
reformatted. The shared **receptor preparation** path has an additional numeric
repair stage that accepts representable exponent/integer fields.

## Raw receptor preparation and model changes

[`prepare_receptor_pdbqt`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/io/receptor_prep.py#L579)
uses a validated preparation wrapper around the canonicalizer. For existing
PDBQT, its sequence is:

1. Repair representable numeric formatting with `normalize_receptor_input`.
   Exponents and integer-like fields can be rewritten if they parse and fit;
   missing coordinates/charges and other fatal lint errors stop preparation.
   Blank occupancy and B-factor fields get defaults of 1.00 and 0.00.
2. Validate the rigid receptor before discarding context: reject multiple
   `MODEL` structures, `ROOT`/`BRANCH`/`TORSDOF`/`BEGIN_RES` records, missing
   atom records, more than 99,999 atoms, unsupported receptor types and blank
   coordinates/charges. Donor hydrogen presence is deferred at this point.
3. Canonicalize to resolve alternate conformations before any recovery.
4. If no `HD` or `HS` atom exists, attempt donor-hydrogen recovery with Open Babel.
5. Canonicalize and validate the final candidate, then replace the destination
   only after preparation succeeds.

These extra checks and recovery are implemented in
[`receptor_checks.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/io/receptor_checks.py).
They are not automatically supplied by the direct `receptor canonicalize`
command or the generic `pdbqt check` command.

Raw PDB is sanitized and converted to PDBQT, then numerically normalized,
validated and canonicalized. MOL2 uses Open Babel by default, retaining its
connectivity for conversion without the PDB sanitizer. A gzip input is unpacked
in temporary storage before processing. Use `--prepare-command` on MolGuard or
`--receptor-prepare-command` on the pipeline for a custom converter.

The PDB sanitizer performs these additional operations:

| Operation | Current policy |
| --- | --- |
| Alternate conformations | Choose per residue: A if present, otherwise the alphabetically first label; retain unlabeled atoms and clear retained labels |
| Blank chain IDs | Assign synthetic chain IDs to blank-chain segments separated by `TER` |
| Occupancy/B-factor fields | Replace missing or unparseable values with 1.00 / 0.00 |
| Element fields | Infer missing elements from atom names and recognized ion residues |
| Solvent and hydrogens | Remove recognized solvent residues and hydrogen records before conversion |
| Residue aliases | ASH/ASQ/DIC/DID→ASP; GLO→GLU; GLZ→GLY; HIZ/HIY→HIS; LEV→LEU; MEU/MSE→MET |
| Selenomethionine | Convert its SE atom to MET SD sulfur |
| Ions | Preserve recognized monoatomic ions and convert their `ATOM` records to `HETATM` |
| Other nonstandard residues | Remove nonstandard `ATOM` residues outside the sanitizer's supported set, with a model-change warning |

For PDB, the default converter is an available Meeko receptor command with
`-i`, `--allow_bad_res` and `--default_altloc A`; Open Babel is the fallback.
For MOL2, the wrapper chooses Open Babel instead of Meeko's PDB reader. The
[Meeko receptor documentation](https://meeko.readthedocs.io/en/develop/cli_rec_prep.html)
describes its preparation interface. Open Babel's `-xr` requests rigid PDBQT
output; `-xn` preserves atom names. The fallback also requests hydrogen addition
and Gasteiger partial charges. Its rigid-output and name options are documented
in the [PDBQT format options](https://openbabel.org/docs/FileFormats/AutoDock_PDBQT_format.html).
A custom converter can be supplied through the preparation options.

If Meeko fails with recognized excess inter-residue bond errors and the command
does not already specify residue deletion, the wrapper retries once with the
affected residue IDs in `--delete_residues`. It prints a prominent warning.
Solvent removal, alias remapping, SE→S replacement, discarded alternate atoms
and residue deletion can change the structural or chemical model. Record and
review those changes separately from numeric formatting.

The conversion subprocess gets a default hash seed and one-thread settings for
common numerical libraries only when those environment variables are unset.
Existing settings are respected. These defaults reduce avoidable variation;
they do not guarantee identical chemistry across converter versions.

## Donor-hydrogen recovery and provenance

Recovery is triggered only when an existing prepared receptor has **no HD/HS
donor hydrogens**. It writes a heavy-atom PDB intermediate, lets Open Babel
perceive bonds, adds hydrogens and recalculates atom types/Gasteiger charges.
The resulting candidate must validate and preserve the original heavy-atom
identity, multiplicity, elements and coordinates at PDB precision.

On successful recovery, the original source bytes are backed up under
`.molguard-backups/<source-name>.<source-sha256>.bak`. A sibling `<stem>.prep.json`
records the method, command, converter diagnostics, hashes, backup path,
preserved heavy-atom count and donor-hydrogen count. Failed recovery does not
replace the source. Already donor-bearing inputs do not trigger recovery.

This geometry-based reconstruction cannot prove bond orders or protonation and
does not reconstruct missing heavy atoms. A partially hydrogenated receptor can
remain incomplete. Review the source chemical states, converter warnings and
model changes; the backup and report support traceability, not chemical approval.

## Where the pipeline uses MolGuard

These operations run automatically in the normal workflow; separate manual
commands are useful for inspection and controlled preparation.

| Stage or entry point | MolGuard operation |
| --- | --- |
| `ultidock setup` / setup stage of `ultidock run` | Discover receptor inputs in `MACRO_MOL_DIR`, prepare/canonicalize with stamp `SETUP`, and scan PDBQTs for lint issues |
| Docking processor initialization | Shared receptor preparation before site/grid loading, with stamp `PIPELINE`; preparation failures stop grid generation |
| Each ligand worker, before receptor/site iteration | Normalize that ligand PDBQT in place; skip a ligand whose normalization fails |
| `benchmarks/prepare_receptor_pdbqt.py` | Use the same receptor-preparation implementation, with default stamp `BENCHMARK` |
| Explicit `molguard receptor prepare` | Shared preparation with default stamp `MOLGUARD` |

Existing same-stem PDBQT is preferred to raw receptor inputs unless preparation
is forced. `--receptor-prep-mode off` disables shared receptor preparation;
skipping setup alone does not disable the docking processor's preparation.
The worker's ligand normalization is independent of that receptor option.

## Reproducibility and manual inspection

For reproducible output bytes, supply a fixed stamp and retain the MolGuard and
converter versions:

```bash
molguard pdbqt check receptor.pdbqt
molguard receptor canonicalize receptor.pdbqt -o receptor_canon.pdbqt --timestamp FIXED
molguard pdbqt normalize ligand.pdbqt -o ligand_normalized.pdbqt
```

The direct receptor canonicalization command defaults to the current UTC time.
Different timestamps change the remark and therefore the hash, even when atom
records match. `SETUP` and `PIPELINE` likewise give different file hashes.
Matching output bytes also require consistent atom labels, ordering of equal
sort keys, formatting and software versions.

SHA-256 identifies the file representation and helps detect changed inputs;
its definition is provided by [NIST FIPS 180-4](https://doi.org/10.6028/NIST.FIPS.180-4).
It is not a chemical identifier or evidence of biological validity. Keep source
structures, prepared files, commands and warnings together. The automatic
backup/report applies to successful donor-hydrogen recovery; ordinary
canonicalization does not create a general source backup.

## Implementation sources and checks

| Source | What it establishes |
| --- | --- |
| [`pdbqt.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/io/pdbqt.py) | Lint rules, ligand preservation, receptor altloc/sort/write/hash algorithm |
| [`fixedfmt.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/io/fixedfmt.py) | Field widths, decimal precision, zero handling and B-factor policy |
| [`receptor_prep.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/io/receptor_prep.py) | Sanitization, converter selection, retry policy and input discovery |
| [`receptor_checks.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/io/receptor_checks.py) | Numeric repair, rigid-receptor validation, donor recovery and heavy-atom preservation |
| [`setup.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/docking/setup.py#L230), [`dock_v02.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/docking/dock_v02.py#L327) | Automatic calls, fixed stamps and worker normalization |
| [`test_pdbqt.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/tests/test_pdbqt.py), [`test_fixedfmt.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/tests/test_fixedfmt.py), [`test_receptor_prep.py`](https://github.com/taka78/ultidock/blob/f617684c4de0363f30c54095ca5dc6cd19a99e7d/molguard/tests/test_receptor_prep.py) | Existing regression checks for formatting, fixed-stamp idempotency, torsion preservation and shared preparation/recovery |

Run those existing checks from an environment with the development dependencies:

```bash
python -m pytest molguard/tests/test_fixedfmt.py molguard/tests/test_pdbqt.py molguard/tests/test_receptor_prep.py
```

For practical input staging, continue to [Automatic receptor preparation](../user-guide/preparing-receptor.md)
or [Ligand inputs](../user-guide/preparing-ligands.md).
