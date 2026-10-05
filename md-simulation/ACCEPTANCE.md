# Native acceptance — 2026-10-05

The soluble software handoff passed after the fixes in this working tree.
A subsequent membrane software test also passed using separately downloaded,
reviewed-for-testing POPC and OPM glycophorin inputs. The repository's SERT
membrane JSON remains a template; it was not treated as a reviewed SERT system.

## Results

| Check | Result |
| --- | --- |
| Fresh network clone of `gmx-beta` | Started at `d0b644f6bae5a263b56e7fab00e6c7d7253ec99d`; fixes below were subsequently installed and tested |
| Clean installation | Python 3.14.4 venv without system site packages; regular, non-editable installation; imports verified outside either checkout |
| `ultidock doctor` | Docking, receptor preparation, GROMACS, ACPYPE, Open Babel and AmberTools reported together |
| Ordinary SERT example | Five successful pocket outputs, 45 poses, 45 matching database rows, 21 filtered analysis rows with binding-site IDs; empty final failure report |
| Real soluble case | T4 lysozyme L99A / benzene, PDB 181L: native Vina docking → selected output pose → GAFF2/AM1-BCC parameterization → build → EM → NVT → NPT → production |
| Deliberately bad ligand | A second benzene candidate with an empty atom map docked successfully, failed MD preparation, and was reported; the valid candidate completed MD |
| Interrupted production | SIGKILL of the runner and native mdrun process group after checkpoint step 80 / 0.16 ps; resumed using `-cpi -append` and unchanged TPR; completed step 2000 / 4 ps |
| Integrity | A 0.25 Å coordinate change in each prepared pose/receptor snapshot was rejected before execution; a changed original docking output was rejected before preparation; original bytes restored |
| Full suite with native smoke enabled | **211 passed, 2 skipped**, 71.34 seconds; skips require unavailable external Holo4K forensic evidence |
| Membrane native acceptance | **Passed for the documented POPC/glycophorin software fixture**; SERT-specific scientific validation remains pending |

The SERT run initially exposed an incompatible system AutoGrid binary and then
missing Vina runtime libraries in the locally extracted test installation.
After correction, the staged example was continued through its ordinary
pipeline. No pre-existing receptor maps or docking results from the developer
checkout were used. The final manifest, pose hashes, database counts, analysis
site IDs and failure CSV were cross-checked.

## Membrane follow-up

The [reproducible membrane test](MEMBRANE_TESTING.md) uses published Slipids
2016 POPC parameters and bilayer coordinates, OPM-oriented 1AFO chains A/B,
and a synthetic benzene pose. This exercises the membrane builder and native
MD execution, not the biological validity of a ligand-binding model. The
default SERT template retains its unset review flags.

The first completed native run retained 92 of 128 lipids, with 46 per leaflet;
added 27 NaCl pairs plus four neutralizing chloride ions; converged EM in
1,085 steps; and completed 2 ps each of NVT, NPT and production. No LINCS
warnings or `-maxwarn` bypass occurred. Native assertions also verify transformed
ligand coordinates, solvent exclusion, salt counts and semi-isotropic coupling.

Review found and fixed periodic pore-water distances and leaflet classification
for wrapped lipid molecules. The builder now rejects a missing retained leaflet,
invalid slab bounds and broken solvent groups. Regression tests cover these
guards, whole-lipid deletion, coordinate/topology mismatches and invalid cells.
Saved inputs and run evidence live under
[`md-simulation/workspace/membrane-acceptance-20261005`](workspace/membrane-acceptance-20261005/),
with the checksum-verified fixture in
[`md-simulation/workspace/membrane-fixture-20261005`](workspace/membrane-fixture-20261005/).

## Fixes found during acceptance

- Rebuild AutoGrid when its advertised `MAX_MAPS` is below the 20 maps required
  by CaV-EMPS. Doctor now reports that limitation. Remove a managed executable
  symlink before rebuilding so its system target is not overwritten.
- Use Meeko's built-in `--read_pdb` reader for automatic receptor preparation;
  the previous `-i` path required undeclared ProDy.
- Isolate selected-ligand chemical-input failures during MD preparation. Save
  JSON/CSV reasons, retain them in the job and execution summary, and mark the
  handoff partial. Valid selected ligands continue; an entirely invalid batch
  returns failure with a report path. Shared receptor/integrity errors remain
  fatal.
- Import SciPy only where membrane construction needs it, so `ultidock md
  --help` also works in a base installation without the optional MD packages.
- Set the mdrun OpenMP environment and command-line budget consistently from
  the protocol. Native testing caught an inherited `OMP_NUM_THREADS=2`
  conflicting with a one-thread protocol.

## Soluble system and limits

The source structure was downloaded from
[RCSB PDB 181L](https://www.rcsb.org/structure/181L). The software test used
chain A's resolved residues 1–162; the unresolved terminal ASN163/LEU164 were
not reconstructed. Crystal waters, HED and chloride were excluded. GROMACS
selected delta-protonated HIS31 and charged termini. The force field was
Amber99SB-ILDN with TIP3P, a dodecahedral box, 0.15 M salt plus neutralization,
300 K, and a 2 fs timestep. The final system contained 38,710 atoms.

Input preparation required retaining the original protein heavy-atom order:
the GROMACS-reordered intermediate caused Meeko to drop serines, which the
MD receptor check rejected. A terminal oxygen was added with pdb2gmx, and
equivalent terminal oxygen labels in the PDBQT were reconciled with the PDB.
All 1,290 protein heavy-atom identities and coordinates then matched before
docking. This was a reviewed preparation of this test input, not a relaxation
of the integrity checks.

EM converged in 711 steps, with Fmax 939.2157 kJ mol⁻¹ nm⁻¹ against the
1,000 threshold. NVT and NPT each ran for 2 ps; production ran for 4 ps.
Recorded energies were finite and the dynamics logs contained no LINCS
warnings. These durations demonstrate software execution and checkpoint
continuation; they do not establish equilibration, binding stability or a
scientifically validated production protocol. The default analysis filters
retained no benzene rows; its 18 raw docking poses remained in the database
and the MD selection correctly used the scored docking outputs independently
of those analysis filters.

For the interruption check, an explicit `--gmx` wrapper forwarded every
command to native GROMACS and added `-cpt 0.01` to mdrun. This made a checkpoint
available during the short test without changing its force field or dynamics
parameters. Recovery before the first checkpoint and graceful signal handling
were not separately exercised by this native interruption check.

## Environment and saved evidence

The isolated installation used Ubuntu 26.04.1, GROMACS 2025.4, ACPYPE 2026.9.4
with its bundled AmberTools executables, Open Babel's wheel, and Meeko 0.8.0.
Native GROMACS, Vina and required libraries were downloaded as distribution
packages and extracted locally rather than installed system-wide. No Conda
environment was used. Standard OS tools and libraries were available through
an explicit PATH; Python user-site packages and the inherited shell environment
were excluded. Every harness command records its environment and exit status.

- [Local logs, commands and verification records](../workspace/acceptance-20261004/)
- [Completed soluble MD job](workspace/acceptance-20261005-soluble/)
- [Full suite log](../workspace/acceptance-20261004/logs/full-suite-final.log)
- [Docking verification](../workspace/acceptance-20261004/docking-verification.json)
- [Checkpoint progression](../workspace/acceptance-20261004/checkpoint-progress.json)
- [Integrity checks](../workspace/acceptance-20261004/integrity-results.json)
- [Molecular assumptions and energy samples](../workspace/acceptance-20261004/soluble-review.json)

These workspace artifacts are intentionally ignored by Git. Their recorded
absolute paths describe the original `/tmp/ultidock-acceptance-20261004` test
installation. The completed MD job was also copied into the repository's MD
workspace; all immutable-input, workflow and completed-output hashes were
verified after copying. This report and the fixes can be committed separately;
the existing remote `gmx-beta` tag has not been moved.
