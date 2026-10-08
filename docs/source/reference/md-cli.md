# MD CLI reference

`ultidock md --help` and `ultidock md COMMAND --help` show the options provided
by `md-simulation/workflow.py`. Examples below use placeholder paths; the
[workflow guide](../user-guide/molecular-dynamics/docking-to-md.md) explains the
inputs and generated job directory.

## Commands

| Command | Purpose |
| --- | --- |
| `ultidock doctor` | Report docking, receptor preparation and MD dependency paths, Python environment and workspaces |
| `ultidock md check --config FILE` | Validate protocol, receptor and native MD tools before docking |
| `ultidock md select --docking-dir DIR --receptor-id ID --engine ENGINE` | Print ranked distinct-ligand poses as JSON |
| `ultidock md prepare --config FILE --docking-dir DIR` | Validate and snapshot a new job; start no native build or dynamics |
| `ultidock md run --config FILE --docking-dir DIR` | Prepare a new job and run through NPT by default |
| `ultidock md run JOB_DIR` | Resume an existing job through NPT by default |

### `md check`

Required: `--config FILE`.

| Option | Effect |
| --- | --- |
| `--receptor-dir DIR` | Check `<receptor_id>.pdbqt` in the actual staged receptor directory |
| `--ligands-dir DIR` | Also validate candidate chemistry/maps against staged PDBQT ligands |
| `--engine {vina,adgpu}` | Verify the protocol matches the selected docking engine |
| `--inputs-only` | Skip native dependency availability checks |
| `--work-dir DIR` | Check that an intended new MD destination is valid and does not already exist |
| `--gmx`, `--acpype`, `--obabel` | Executable overrides; defaults are those command names |

### `md select` and `md prepare`

`select` requires `--docking-dir`, `--receptor-id` and `--engine {vina,adgpu}`.
`--top N` defaults to five distinct ligand IDs and must be positive.
`prepare` requires `--config` and `--docking-dir`; its count comes from protocol
`top`, and `--work-dir` selects a new destination under `md-simulation/`.

Both accept `--legacy-single-receptor` for historical filenames lacking a
receptor ID. Use it only after confirming that the directory contains results
for exactly the declared receptor. Neither command accepts `--pose-manifest`;
use the direct `md run --config` route to restrict preparation to a recorded run.

### `md run`

Use either `--config FILE --docking-dir DIR` for a new job or positional
`JOB_DIR` to resume; these forms cannot be combined.

| Option | Effect |
| --- | --- |
| `--through {prepare,build,em,nvt,npt,production}` | Last stage; default `npt`; `prepare` applies only to a new job, production requires resuming a reviewed job |
| `--equilibration-reviewed` | Required when resuming through production |
| `--gmx`, `--acpype`, `--obabel` | Executable overrides for new or resumed execution |
| `--work-dir DIR` | New-job destination under the active `md-simulation/`; use positional `JOB_DIR` for an existing one |
| `--pose-manifest FILE` | New-job selection limited to the listed docking invocation |
| `--receptor-dir DIR` | Bind a new job to the receptor staged for that docking invocation |
| `--result-file FILE` | Write new-job handoff status and job path as JSON |
| `--legacy-single-receptor` | Permit reviewed historical single-receptor output names for a new job |

The manifest, receptor directory, result file and legacy options are for the
new-job form. Resume uses the saved inputs/provenance. Shared integrity and
dependency failures stop before simulation; per-system failures are retained
in `run_summary.json` and return a failure exit status after processing.

## Continue from public docking commands

`run`, `cavity`, `blind`, `known-site`, `fpocket`, `p2rank`, `run fpocket/p2rank`
and the SERT/4COF example runners accept:

| Option | Effect |
| --- | --- |
| `--md-config FILE`, alias `--md FILE` | Enable continuation using a completed JSON protocol |
| `--md-through {prepare,build,em,nvt,npt}` | Last stage for the new job; default `npt` |
| `--md-work-dir DIR` | New MD directory under the active `md-simulation/` |
| `--md-gmx`, `--md-acpype`, `--md-obabel` | MD executable overrides |

Other `--md-*` options require `--md-config`. The docking backend must match
protocol `engine`; explicit GPU mode pairs with `adgpu`, CPU with `vina`.
`--skip-analysis` does not disable MD. Use a positional MD job directory to
continue into production after reviewing equilibration.
