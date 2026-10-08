# MD troubleshooting

Start with `ultidock doctor`, the docking handoff `*.md.json`, and the first
failed `*.command.log` in the printed MD job. Preserve the job and its source
docking artifacts while diagnosing it. See [MD installation](../getting-started/md-installation.md)
and [restart behavior](../user-guide/molecular-dynamics/running-and-results.md).

## Dependencies and handoff

| Symptom | Check and action |
| --- | --- |
| `MD dependencies: incomplete` | Activate the same Python environment, source AmberTools `amber.sh` and any GROMACS `GMXRC`, then check every tool path in doctor |
| `antechamber`, `parmchk2`, `tleap` or `sqm` missing | Install/load AmberTools; a platform-specific ACPYPE wheel does not guarantee that Ultidock finds these commands on `PATH` |
| ACPYPE cannot import Open Babel | Install its Python bindings in ACPYPE's environment, such as `openbabel-wheel` for the pip route |
| Native tool fails despite `[OK]` in doctor | Availability checks do not execute the full toolchain; inspect its command log for libraries, version or parameter errors |
| Protocol `engine` does not match docking | Use `adgpu` for CUDA/OpenCL or `vina` for CPU, with matching docking `--mode` |
| Handoff says `skipped` | Read its reason: MD preflight may have failed, all docking cases may have failed, or the requested receptor/engine may have no successful outputs |
| CSV is empty but poses exist | MD selection reads scored engine outputs, independently of the analysis filters; check the manifest and requested receptor/engine |
| Missing AutoDock-GPU coordinates | Retain the paired XML and DLG; a `-best.pdbqt` alone cannot identify the exact scored run |
| `XML/DLG energy mismatch` or missing model | Restore intact paired outputs from that run or rerun docking; do not substitute another ligand/run's coordinates |

## Molecular inputs

| Symptom | Check and action |
| --- | --- |
| `receptor_reviewed` is false | Review completeness, chemical states, termini, disulfides and required components, then complete the protocol |
| Receptor heavy atoms do not match | Compare the protein PDB with the exact docking PDBQT; restore the shared frame/identity or prepare and dock the revised structure |
| Selected poses used another receptor preparation | Use a manifest/directory belonging to that exact receptor; matching filenames alone are insufficient |
| Ligand missing from `ligands` | Key entries by selected ligand filename stem and supply its SDF, integer charge and complete atom map |
| Atom map is incomplete or nonunique | Map every PDBQT heavy-atom serial to its original zero-based SDF atom index once; do not guess from atom names or proximity |
| Formal charge, stereochemistry or bond lengths disagree | Recover the chemical state actually docked from the original preparation; check the SDF/map and input coordinates |
| One ligand is skipped during preparation | Read `preparation_failures.json`/`.csv`; other valid selected ligands continue, and lower-ranked candidates are not substituted |
| No selected ligands pass preparation | The command fails with the saved report path; correct chemistry inputs and prepare a new job |
| Unsupported force-field/water combination | Use the supported Amber/GAFF2 combinations and required external `.ff`/water templates listed in [Protocol reference](md-protocol.md) |
| Missing/renamed ligand topology atoms or incompatible defaults | Inspect AmberTools/ACPYPE logs, atom ordering, charge and protein/lipid parameter compatibility; correct inputs and create a fresh job |

## Membrane systems

| Symptom | Check and action |
| --- | --- |
| Orientation review or compatibility missing | Review the actual seed/parameters, common protein–ligand transform and protein force-field match |
| Lipid coordinate/topology mismatch | Check ordered molecule blocks, counts, atom names and atoms per molecule against the lipid-only seed |
| Nested lipid includes rejected | Supply flattened include files with all parameter dependencies |
| Tilted cell or dodecahedron rejected | Use an orthorhombic periodic seed with membrane normal Z |
| Missing retained leaflet | Inspect insertion, periodic molecule positions and clash threshold; both leaflets must retain lipids |
| Incorrect slab/pore hydration | Slab and pore values are in nm in the seed frame; the rigid protein/ligand transform acts on Å coordinates |
| Cutoff exceeds padding or switching is invalid | Use the lipid authors' nonbonded settings and adequate periodic clearance; switching needs `0 < rvdw-switch < rvdw` |

## Execution and continuation

| Symptom | Check and action |
| --- | --- |
| Destination already exists | Resume with `ultidock md run JOB_DIR`, or choose a new preparation directory |
| Work directory outside `md-simulation/` | Use a destination under the MD directory in the active workspace printed by doctor |
| Production refused | Run through NPT, inspect equilibration, then resume the saved job with `--through production --equilibration-reviewed` |
| GROMACS warnings stop `grompp` | Inspect and correct the actual topology/geometry/settings; the runner does not pass `-maxwarn` |
| EM did not converge | Inspect `em.log`, clashes and parameters; a finished process alone is not accepted as converged minimization |
| System build was interrupted or failed | Correct the inputs and prepare a fresh job; topology/solvent edits are not guessed at during recovery |
| Dynamics interrupted | Repeat the saved endpoint command; the runner reuses a matching TPR/checkpoint with `-cpi -append` |
| TPR, input or completed output changed | Restore the immutable original or prepare a new job; do not edit saved hashes to bypass validation |
| Workflow code changed after preparing the job | Prepare a fresh job with the current workflow; its code hashes are part of provenance |
| `This MD job is already running` | A job lock prevents simultaneous execution of the same job; inspect the current process before resuming |
| One complex fails and others finish | Inspect `run_summary.json` and that complex's command/state logs; the batch returns failure after retaining successes |

For a bug report, include the branch/commit, executable versions, protocol,
handoff status, failed command log and `run_summary.json`, with shareable
minimal inputs. State whether failure happened in selection, chemistry,
system build, minimization, dynamics or restart.
