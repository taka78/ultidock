# CLI and script reference

`ultidock --help` lists commands. A pip install exposes `ultidock` and
`molguard` on `PATH`; from an uninstalled checkout, use
`python -m cli.ultidock` and `python -m cli.molguard` with matching Python
dependencies. The [installation guide](../getting-started/installation.md)
also documents the Ubuntu `/usr/bin/python3` path.

## Workflow commands

| Command | Purpose |
| --- | --- |
| `ultidock doctor` | Show active workspace and docking/MD tools, including AutoGrid map-capacity diagnostics |
| `ultidock setup [options]` | Create directories, detect/build tools and generate `docking/config.py` |
| `ultidock run [options]` | Setup, extraction, docking and analysis; optional `--md-config PROTOCOL` continues into MD |
| `ultidock md COMMAND` | Validate, select poses, prepare, run or resume GROMACS jobs; see [MD CLI](md-cli.md) |
| `ultidock known-site --center x,y,z [options]` | Dock in a supplied expert box |
| `ultidock cavity [options]` | Propose receptor-only CaV-EMPS boxes, then dock |
| `ultidock blind [options]` | Whole-receptor box baseline |
| `ultidock fpocket` / `ultidock p2rank` | Predict local pockets and dock |
| `ultidock run fpocket` / `ultidock run p2rank` | Equivalent pocket forms; choose the sole PDBQT in the default receptor directory unless `--receptor PATH` is passed |
| `ultidock pocket-box --method METHOD --receptor FILE --output TSV` | Create fpocket/P2Rank boxes without docking |
| `ultidock profile-receptors` | Generate receptor geometry `.config.toml` sidecars |
| `ultidock report DIRECTORY` | Build Markdown/HTML and visualization helpers |
| `ultidock clean [options]` | Remove selected generated artifacts; inspect `ultidock clean -- --help` before deletion |

`known-site` is a manual-box baseline; `cavity` generates receptor-only
site hypotheses; `blind` uses one broad region. Site IDs are output names,
not rankings. Pocket methods install only the chosen local tool if absent.
For `fpocket` and `p2rank`, `--tool PATH` overrides the local executable;
`--autosites` defaults to 6, `--box-size` to 35 Å, and
`--grid-spacing` to 0.375 Å. The box center is the predictor's reported
center (fpocket uses its pocket-coordinate centroid); no extra pocket
padding is added. `sites.tsv` and raw predictor output stay with the run.

## Examples and benchmarks

| Command | Purpose |
| --- | --- |
| `ultidock example list` | Discover bundled examples |
| `ultidock example run NAME [arguments]` | Stage a declared example in a new workspace |
| `ultidock benchmark download-dude` | Fetch DUD-E target files |
| `ultidock benchmark cavity-recovery` | Evaluate receptor-only site localization |
| `ultidock benchmark site-prediction` | Run COACH420/HOLO4K predictor comparisons |
| `ultidock benchmark dude-docking` / `full-dude` | Run active/decoy docking evaluation |
| `ultidock benchmark plots` | Plot supported benchmark summaries |

See [Bundled examples](../tutorials/examples.md) and
[Benchmarking](../benchmarks/index.md) for complete commands and outputs.

## Low-level checkout scripts

| Script | Role |
| --- | --- |
| `python3 docking/run.py [options]` | Primary setup → extract → dock → analyze entry point |
| `python3 docking/setup.py [options]` | Setup only; most flags shared with `run.py` |
| `python3 docking/dock_v02.py` | Internal receptor/site/ligand docking stage |
| `python3 docking/extract.py` | Archive extraction and `vina_split` staging |
| `python3 docking/clean.py -y --all` | Explicit full reset of generated files |
| `python3 docking/profile_receptors.py` | Geometry sidecar generation |
| `python3 docking/analyse_docking_results.py` | Analysis from SQLite |
| `python3 benchmarks/download_dude.py` | DUD-E downloader |
| `python3 benchmarks/cavity_recovery_benchmark.py` | Site-recovery evaluation with target-level `--jobs` |

Use the public CLI for normal work. A clean `--all` request removes
compiled tools, maps, downloaded ligands, configuration and results; inspect
its scope first. `ultidock clean -y` alone is limited to build files and
caches. See [Advanced setup and pipeline controls](../getting-started/setup-and-run.md).

## Pipeline options and help forwarding

`ultidock run -- --help` forwards help to the underlying script, while
`ultidock run --help` displays wrapper help. `python docking/run.py --help`
shows the same underlying options. Wrapper commands forward ordinary pipeline
flags after their own options. Key flags are `--mode`, `--skip-wget`,
`--skip-setup`, input/output directory overrides, receptor preparation
controls, grid settings and the Vina sampling controls. The full list and
defaults are in [Configuration](configuration.md).

SERT accepts an optional positional `p2rank` or `fpocket` method:

```bash
ultidock example run quickstart -- --help
ultidock example run quickstart --dry-run
ultidock example run sert-escitalopram -- --help
ultidock example run sert-escitalopram p2rank --mode gpu
```

The quickstart dry run only prints planned demonstration outputs. Native SERT
and 4COF example dry runs stage inputs and print the pipeline command.
High-level `known-site` and `cavity` previews can write a run directory and
`run_config.yaml`; pocket-method dry runs can execute prediction. Read each
command's help rather than assuming every preview has identical behavior.

## Molecular dynamics

Add `--md-config PROTOCOL` (alias `--md PROTOCOL`) to a docking pipeline to
continue using its current successful-output manifest. The default endpoint is
NPT; `--md-through` can stop at prepare, build, EM, NVT or NPT. Use
`--md-work-dir` for a new job directory under the active `md-simulation/` tree.
For existing results, `ultidock md run --config PROTOCOL --docking-dir PATH`
selects, prepares and simulates. Existing jobs resume with `ultidock md run JOB`.

See [MD installation](../getting-started/md-installation.md),
[MD CLI](md-cli.md), [protocol fields](md-protocol.md) and
[docking-to-MD walkthrough](../user-guide/molecular-dynamics/docking-to-md.md)
for the required inputs and full commands. Production additionally requires
`--equilibration-reviewed` after inspection of the system and equilibration.

## MolGuard commands

| Command | Result |
| --- | --- |
| `molguard doctor` | Version and optional converter diagnostics |
| `molguard pdbqt check FILE` | Receptor/ligand PDBQT lint, including fixed-column, numeric and atom-type errors |
| `molguard pdbqt normalize FILE -o OUT` | Reformat ligand numeric fields without changing its torsion tree |
| `molguard receptor canonicalize FILE -o OUT` | Deterministically sort/renumber receptor PDBQT and report a SHA-256 digest |
| `molguard receptor prepare FILE -o OUT` | Canonicalize PDBQT without adding atoms, or sanitize/convert a raw receptor |
| `molguard grids check MAPS.fld` | Check missing/all-zero/nonfinite maps and type mismatches |

Receptor preparation and ligand normalization have different chemistry and
record handling. Use the command for the molecule's actual role. The
[receptor guide](../user-guide/preparing-receptor.md) describes the conversions,
validation and preservation of source inputs.
