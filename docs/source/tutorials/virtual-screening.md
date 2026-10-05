# Virtual screening

A screen extends a validated single-target workflow to a traceable ligand library.
First run a small pilot containing representative ligand sizes, flexible molecules
and any unusual atom types expected in the full library.

## 1. Prepare the library

Assign unique filenames to chemical states and retain a table mapping each filename
to source compound, stereochemistry, protonation and preparation settings. Validate
the resulting PDBQTs with MolGuard. Keep only intended inputs in the ligand directory.

## 2. Fix the protocol

Choose a receptor, site definition, engine and sampling settings. For a known site,
substitute your actual coordinates and paths in this preview:

```bash
ultidock known-site --center=12.3,4.5,-6.7 --box-size 25 \
  --mode cpu --vina-cpu 2 --vina-exhaustiveness 8 --vina-seed 42 \
  --skip-wget --macro-mol-dir /absolute/screen/receptors \
  --ligands-dir /absolute/screen/ligands \
  --output-dir /absolute/screen/pilot --dry-run
```

The coordinates are illustrative, not a D2 binding-site claim. Repeat without
`--dry-run` after reviewing your own box. A seed controls one sampling setting;
keep versions and all other parameters with the run as well.

## 3. Scale after checking the pilot

Verify expected ligand counts, job completion, raw pose files and score associations.
Estimate storage and runtime from the pilot. Use separate ligand shards and job
workspaces for [HPC screening](hpc-screening.md), with identical preparation and
search parameters when comparing shards.

## 4. Review candidates

Merge exports by receptor, site, ligand state and engine metadata; do not join on
an ambiguous compound name alone. Inspect poses, retain failures, and evaluate
ranking with independently established labels if available. A top-ranked docking
score alone is not an experimentally validated hit.
