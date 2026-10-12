# How MolGuard canonicalizes molecular inputs

MolGuard gives Ultidock a consistent representation of prepared docking inputs.
For a **receptor**, canonicalization selects alternate-location records, orders
atoms, replaces serial numbers and writes a standardized PDBQT. For a **ligand**,
normalization rewrites numeric fields while preserving its atom order and torsion
records. Raw receptor conversion is a separate preparation step.

| Operation | Input | Main effect |
| --- | --- | --- |
| `molguard pdbqt check` | Receptor or ligand PDBQT | Report formatting errors and warnings without writing |
| `molguard receptor canonicalize` | Prepared rigid-receptor PDBQT | Select altloc, sort/renumber atoms, reformat and return an output hash |
| `molguard pdbqt normalize` | Prepared ligand PDBQT | Reformat atom records while preserving atom serials/order and other records |
| `molguard receptor prepare` | PDB, MOL2 or PDBQT, including gzip inputs | Sanitize PDB when applicable, invoke a converter if needed, then canonicalize |

This page describes the implementation on `gmx-dev` at commit
`7862e2cd1515e278c8235050aefa92220f213fdc`. The implementation links below
identify the policies used by this version.

## Why these steps are needed

PDBQT extends PDB coordinates with partial charges, AutoDock atom types and
ligand torsion records. Its numeric layout and torsion-tree conventions are
specified in the [AutoDock 4.2.6 User Guide, Appendix I](https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf).
The [wwPDB coordinate specification](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)
defines the underlying coordinate fields, alternate-location identifiers and
insertion codes.

The bundled [AutoGrid parser](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/docking/AUTODOCK_GPU_DIR/autogrid/mainpost1.28.cpp#L875)
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
[`canonicalize_receptor`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/io/pdbqt.py#L457).

### 1. Check the input

`pdbqt_check` decodes the file as ASCII and checks the numeric slices of each
`ATOM`/`HETATM` record. It collects truncated fields, exponent notation, missing
decimal points, unparseable numbers, NaN/Inf and missing atom-type fields.
Short atom records and B-factors above 99.99 produce warnings. A lint error
stops canonicalization or normalization before the output is written.

This checker covers a limited format contract. It does not validate bond orders,
protonation, structural completeness, torsion-tree consistency or membership of
each nonblank atom type in an allowed parameter set. Blank numeric fields can
also fail later when parsed. Receptor atom-type compatibility is checked
separately by [grid preparation](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/docking/make_grids.py#L288).

### 2. Parse atom records and select alternate locations

The canonicalizer parses `ATOM` and `HETATM` records into atom objects. If
nonblank altloc labels exist, it keeps the alphabetically first label found
**across the file**, together with all blank-altloc atoms. Other labeled atoms
are discarded with a warning; retained altloc labels are cleared.

For example, a prepared PDBQT containing labels A and B retains A. This is a
deterministic selection rule, not a highest-occupancy calculation or a review of
which conformation is biologically appropriate. Select the intended conformer
before docking when this choice matters. Raw PDB sanitization uses a different,
per-residue rule described below.

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
[`molguard.io.fixedfmt`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/io/fixedfmt.py).
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

[`pdbqt_normalize`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/io/pdbqt.py#L295)
applies the same lint gate and numeric formatter, but writes each record back in
its original order. Atom serials and nonnumeric atom fields are preserved.
`ROOT`, `ENDROOT`, `BRANCH`, `ENDBRANCH`, `TORSDOF` and other non-atom records
are passed through unchanged apart from line-ending normalization.

Ligand branch references depend on atom serials and record placement, which is
why this path avoids sorting or renumbering. A formatting error such as
`EXPONENT` or `NO_DECIMAL` is **rejected**, not automatically repaired by these
commands. Correct the source writer or regenerate a valid prepared file, then
normalize it. Warning-only B-factor issues can be reformatted.

## Raw receptor preparation and model changes

[`prepare_receptor_pdbqt`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/io/receptor_prep.py#L565)
canonicalizes existing PDBQT without adding atoms or assigning new charges.
Raw PDB is first sanitized, converted to PDBQT and then canonicalized. MOL2
goes through an appropriate converter and canonicalization without the PDB
sanitizer. A gzip input is unpacked in temporary storage before processing.
For MOL2, supply a converter that reads that format through `--prepare-command`
on MolGuard or `--receptor-prepare-command` on the pipeline; the auto-detected
Meeko command uses its PDB reader.

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

The default converter is an available Meeko receptor command with `--read_pdb`,
`--allow_bad_res` and `--default_altloc A`; Open Babel is the fallback. The
[Meeko receptor documentation](https://meeko.readthedocs.io/en/develop/cli_rec_prep.html)
describes its preparation interface. Open Babel's `-xr` requests rigid PDBQT
output, as documented in its [PDBQT format options](https://openbabel.org/docs/FileFormats/AutoDock_PDBQT_format.html).
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

## Where the pipeline uses MolGuard

These operations run automatically in the normal workflow; separate manual
commands are useful for inspection and controlled preparation.

| Stage or entry point | MolGuard operation |
| --- | --- |
| `ultidock setup` / setup stage of `ultidock run` | Discover receptor inputs in `MACRO_MOL_DIR`, prepare/canonicalize with stamp `SETUP`, and scan PDBQTs for lint issues |
| Docking processor initialization | Shared receptor preparation before site/grid loading, with stamp `PIPELINE`; collect receptor-preparation failures |
| Each ligand worker, before receptor/site iteration | Normalize that ligand PDBQT in place; report failed normalization and continue other ligand cases |
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
It is not a chemical identifier or evidence of biological validity. Docking
[pose provenance](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/docking/pose_provenance.py)
and the MD handoff also hash actual receptor/ligand artifacts to check their
identity. Keep source structures, prepared files, commands and warnings together;
automatic in-place preparation does not create a source backup here.

## Implementation sources and checks

| Source | What it establishes |
| --- | --- |
| [`pdbqt.py`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/io/pdbqt.py) | Lint rules, ligand preservation, receptor altloc/sort/write/hash algorithm |
| [`fixedfmt.py`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/io/fixedfmt.py) | Field widths, decimal precision, zero handling and B-factor policy |
| [`receptor_prep.py`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/io/receptor_prep.py) | Sanitization, converter selection, retry policy and input discovery |
| [`setup.py`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/docking/setup.py#L230), [`dock_v02.py`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/docking/dock_v02.py#L354) | Automatic calls, fixed stamps and worker normalization |
| [`test_pdbqt.py`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/tests/test_pdbqt.py), [`test_fixedfmt.py`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/tests/test_fixedfmt.py), [`test_receptor_prep.py`](https://github.com/taka78/ultidock/blob/7862e2cd1515e278c8235050aefa92220f213fdc/molguard/tests/test_receptor_prep.py) | Existing regression checks for formatting, fixed-stamp idempotency, torsion preservation and shared preparation |

Run those existing checks from an environment with the development dependencies:

```bash
python -m pytest molguard/tests/test_fixedfmt.py molguard/tests/test_pdbqt.py molguard/tests/test_receptor_prep.py
```

For practical input staging, continue to [Automatic receptor preparation](../user-guide/preparing-receptor.md)
or [Ligand inputs](../user-guide/preparing-ligands.md). The
[MD chemistry guide](../user-guide/molecular-dynamics/soluble-systems.md)
explains the additional chemical graph and coordinate checks required for MD.
