# Run stages, results and restart

New MD jobs default to restrained NPT. Use the printed job directory to inspect
or resume a job; its source inputs and workflow hashes are checked before
execution. [Docking to MD](docking-to-md.md) explains how that job is created.

## Simulation stages

| Endpoint | What has completed |
| --- | --- |
| `prepare` | Selected poses, chemistry, immutable input snapshots and generated MDP files; no system build or dynamics |
| `build` | Protein and ligand parameters, assembled/solvated system, membrane insertion when requested, and ions |
| `em` | Energy minimization; the runner checks that GROMACS reports convergence |
| `nvt` | Restrained constant-volume equilibration; initial velocities generated once |
| `npt` | Restrained constant-pressure equilibration; continues the NVT checkpoint |
| `production` | Unrestrained dynamics; continues the NPT checkpoint after review |

The baseline constrains bonds involving hydrogen, a timestep at most 2 fs, a single
v-rescale thermostat bath and C-rescale pressure coupling. Pressure coupling is
isotropic for soluble systems and semi-isotropic for membranes. Protein/ligand
heavy-atom restraints are active during NVT/NPT and removed for production.
See [MD methodology](../../scientific-background/molecular-dynamics.md) for what
these choices mean.

## Review before production

Inspect the built topology and coordinates, temperature, density, pressure
trends, ligand pose and contacts. For a membrane, also review packing, leaflet
balance and pore hydration. The configured duration alone does not establish
equilibration. When ready, continue the saved job:

```bash
ultidock md run /path/to/md-job --through production --equilibration-reviewed
```

The review is recorded in `equilibration_review.json` with the NPT checkpoint
hash. Production cannot be requested as the first direct `md run --config`
invocation; prepare and run through NPT, then resume the job directory.

## Find the artifacts

Default jobs are under the active workspace's `md-simulation/workspace/`.
Custom work directories must also stay under its `md-simulation/` directory.
The CLI prints the job path and a command with any executable overrides.

| Location or file | Contents |
| --- | --- |
| `job.json`, `job.sha256.json` | Selected ligand systems, protocol, pose provenance and immutable-input/workflow hashes |
| `submitted_protocol.json`, `inputs/` | Original protocol and snapshotted receptor/force-field/bilayer inputs |
| `docking_manifest.json` | Current-run manifest snapshot when selection was restricted to a docking invocation |
| `preparation_failures.json`, `.csv` | Ligands rejected during MD preparation, with score/site and reason |
| `01_<ligand>/`, `02_<ligand>/`, … | One directory per prepared complex |
| `selected_pose.pdbqt`, `chemical_source.sdf`, `posed_ligand.mol` | Selected docking coordinates and chemically reconstructed ligand |
| `topol.top`, `ligand.itp`, `ligand_atomtypes.itp`, `system.gro` | Built system topology, parameters and coordinates |
| `salt.json`, `membrane_composition.json` | Ion bookkeeping and membrane composition where applicable |
| `*.mdp`, `*.processed.mdp`, `*.tpr` | Stage settings, processed settings and binary simulation input |
| `*.log`, `*.edr`, `*.gro`, `*.cpt`, `*.xtc` | Native logs, energies, final coordinates, checkpoints and sampled trajectories |
| `*_energy.xvg` | Extracted temperature, pressure, density and potential-energy traces |
| `commands.jsonl`, `*.command.log`, `state.json` | Command arguments, native output and stage status per complex |
| `run_summary.json` | Completed/failed systems and retained preparation failures |

MD trajectories and stage state stay in these job files. The docking SQLite
database remains the docking-score store; it is not a trajectory database.
Keep the MD job and docking manifest/output containers together for provenance.

## Resume an interrupted dynamics stage

Repeat the saved job command for the desired endpoint:

```bash
ultidock md run /path/to/md-job --through npt
```

A dynamics checkpoint is resumed using its original TPR with `-cpi` and
`-append`. Completed-stage outputs are checked before reuse. Missing or changed
TPRs, completed outputs, immutable inputs or workflow code stop continuation.
Correct changed inputs and prepare a fresh job instead of editing `job.json`
or generated MDP files in place.

A failed or interrupted **system build** requires a new job because its commands
modify topology and solvent in place. A failed complex stops its stages while
the runner continues other prepared complexes, writes `run_summary.json` and
returns failure when any execution failed. Shared input-integrity or dependency
errors stop before simulation. See [MD troubleshooting](../../reference/md-troubleshooting.md).

## Hardware and performance planning

The current runner processes ligand complexes **sequentially**. Its protocol
`threads` value is passed to `gmx mdrun -nt THREADS -ntomp THREADS` and to
`OMP_NUM_THREADS`. Docking's `ULTIDOCK_WORKERS` and `GPU_SLOTS_PER_DEV` do not
set MD concurrency. The runner uses thread-MPI `gmx` directly; it has no MPI
launcher or automatic distribution of independent MD systems across GPUs.

| Change | Expected effect and useful check |
| --- | --- |
| Increase protocol `threads` | Changes the CPU budget for one simulation; measure throughput within your allocated cores |
| Use a GPU-enabled GROMACS | Compatible tasks may be offloaded by GROMACS; confirm device/task assignment in its log |
| Add solvent, lipids or a larger box | Increases atom count, memory, time and trajectory size |
| Increase `top` or stage durations | Adds separate complexes or timesteps to execute |
| Use faster storage | May help input/output overhead; retain space for trajectories, checkpoints and raw docking artifacts |

GPU support and task assignment depend on the GROMACS build and simulation
settings. Check [GROMACS performance guidance](https://manual.gromacs.org/current/user-guide/mdrun-performance.html)
before choosing hardware. Multiple docking GPUs do not guarantee the same
scaling for a single MD simulation.

Measure **ns/day** from a representative GROMACS run with the same atom count,
settings and hardware. Estimated production days are total planned simulation
nanoseconds divided by measured ns/day. For example, five 100 ns systems at
20 ns/day suggest 25 days of production, plus preparation and equilibration.
This is an estimate; the repository's short native tests are software checks,
not a hardware throughput benchmark.
