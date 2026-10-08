# Docking to MD, step by step

This guide starts with your own docking inputs and a reviewed MD protocol.
Install the [MD dependencies](../../getting-started/md-installation.md) first.
The [docking guide](../start-docking.md) explains `MACRO_MOL_DIR`, `LIGANDS_DIR`
and the optional `ligands.wget` download list.

## 1. Fill in a protocol for your system

Copy the soluble or membrane JSON template from the active workspace's
`md-simulation/examples/`. For example, from an editable checkout root:

```bash
cp md-simulation/examples/soluble.json protocol.json
```

Edit its paths and settings for the actual input. Paths inside the JSON are
relative to **that JSON file**, so moving the template requires updating its
relative paths. The template starts with `receptor_reviewed=false` and an empty
atom map so an unreviewed example cannot start simulation.

Set `receptor_id` to the docking receptor's filename stem. Provide
`receptor_pdb` for a reviewed, protein-only structure and
`docking_receptor_pdbqt` for the exact docking preparation in the same frame.
Fill in the ligand SDF, integer `net_charge` and complete heavy-atom `atom_map`
under each ligand's filename stem. [Soluble systems](soluble-systems.md) explains
these inputs and the supported force fields. For a transporter or other
membrane protein, follow [Membrane systems](membrane-systems.md).

Choose an engine consistent with the docking backend:

| Protocol `engine` | Docking backend | Example |
| --- | --- | --- |
| `adgpu` | AutoDock-GPU on CUDA or OpenCL | `ultidock cavity --mode gpu --md-config ./protocol.json` |
| `vina` | CPU Vina | `ultidock cavity --mode cpu --md-config ./protocol.json` |

Ordinary docking defaults to automatic GPU/CPU selection. An integrated MD
protocol declares one engine, so an explicit matching mode makes the handoff
repeatable. `--mode gpu` requires a detected GPU; multiple detected NVIDIA GPUs
are used by the docking scheduler. Add `--skip-wget` when using only local ligands.

## 2. Check the inputs and tools

```bash
ultidock doctor
ultidock md check --config ./protocol.json
```

`md check` validates the protocol, receptor correspondence, stage settings and
MD tools. To also check candidate ligand chemistry before docking, pass the
actual ligand directory:

```bash
ultidock md check --config ./protocol.json --ligands-dir ./docking/LIGANDS_DIR
```

Candidate-specific chemical failures produce warnings. Supply chemical inputs
for every candidate you want eligible for MD; if a selected ligand lacks a
valid SDF or map, it is skipped and is not replaced by a lower-ranked candidate.
Use `--inputs-only` to check protocol inputs without requiring native MD tools.

## 3. Dock and continue through equilibration

For a selected `ligands.wget` library and an `adgpu` protocol:

```bash
ultidock cavity --mode gpu --md-config ./protocol.json
```

For an already staged local library:

```bash
ultidock cavity --mode gpu --skip-wget --md-config ./protocol.json
```

The same MD options work with `ultidock run`, `blind`, `known-site`, `fpocket`,
`p2rank` and `run fpocket/p2rank`. Ultidock checks MD inputs before docking,
records the current run's successful outputs, selects distinct ligands, prepares
the MD job, and runs through NPT. `--skip-analysis` skips the docking CSV export
but still allows MD. Use `--md-through prepare`, `build`, `em` or `nvt` to stop
earlier. Default job directories are under `md-simulation/workspace/` in the
active workspace; the CLI prints the exact path and resume command.

The `docking-run-*.json` manifest restricts integrated selection to the current
docking invocation. Its sibling `docking-run-*.md.json` records the handoff
status and MD job path. Earlier poses in a reused docking directory are excluded
from this integrated handoff.

## 4. Continue from existing docking results

For a completed screen, use its actual docking output directory and manifest:

```bash
ultidock md run --config ./protocol.json \
  --docking-dir /path/to/run/docking \
  --pose-manifest /path/to/run/docking/docking-run-TIMESTAMP.json
```

Without `--pose-manifest`, this command considers all matching results in that
directory. Restrict a reused directory to one recorded invocation when needed.
The CLI prints the new job directory. A custom `--work-dir` must be a new
directory under the active workspace's `md-simulation/`.

For an inspectable selection/preparation sequence:

```bash
ultidock md select --docking-dir /path/to/run/docking \
  --receptor-id protein --engine adgpu
ultidock md prepare --docking-dir /path/to/run/docking --config ./protocol.json
ultidock md run /path/printed/by/prepare --through npt
```

`select` and `prepare` scan their directory; `--pose-manifest` is available on
the direct `md run --config` route. Selection reads scores and exact coordinates
from engine outputs rather than relying on filtered CSVs or historical database
paths. AutoDock-GPU needs its paired XML and DLG; a `-best.pdbqt` alone is
insufficient. Vina uses the scored `MODEL` in its output PDBQT.

## 5. Review and start production

Inspect the built complex, parameters and NVT/NPT outputs using
[Run stages, results and restart](running-and-results.md). When the system is
ready for production, resume the **existing job directory**:

```bash
ultidock md run /path/to/md-job --through production --equilibration-reviewed
```

The option records your review; it does not calculate whether equilibration
was sufficient. Production duration comes from the snapshotted protocol.
Reuse any `--gmx`, `--acpype` and `--obabel` overrides printed by the first run.
Input or workflow changes require preparing a fresh job.

## If part of the screen fails

Failed docking cases are listed in `docking-run-*.failures.csv`, while successful
cases continue. If all cases fail, analysis and MD are skipped and the pipeline
returns failure. If the requested MD receptor has no successful outputs, the
handoff is marked `skipped` while other receptor results remain available.
An invalid MD preflight disables the handoff while docking continues.

Individual ligand preparation failures are recorded in
`preparation_failures.json` and `.csv`; valid selected ligands continue. An MD
execution failure preserves the job and returns failure after the prepared
systems have been processed. See [MD troubleshooting](../../reference/md-troubleshooting.md)
for the relevant logs and restart rules.
