# CLI reference

`ultidock --help` lists commands; append `--help` to a command for its own options.
From an uninstalled checkout, substitute `python -m cli.ultidock`. MolGuard's
checkout entry point is `python -m cli.molguard`.

| Command | Purpose |
| --- | --- |
| `ultidock doctor` | Report tools, paths and active workspace |
| `ultidock setup [options]` | Prepare tools and generate pipeline configuration |
| `ultidock run [options]` | Run the full pipeline |
| `ultidock known-site --center X,Y,Z` | Dock inside a supplied box |
| `ultidock cavity --autosites N` | Propose CaV-EMPS sites and dock |
| `ultidock blind` | Use a whole-receptor search region |
| `ultidock fpocket` / `ultidock p2rank` | Predict pockets locally and dock |
| `ultidock pocket-box --method METHOD --receptor FILE --output TSV` | Generate boxes without docking |
| `ultidock profile-receptors` | Generate receptor geometry sidecars |
| `ultidock example list` | Discover bundled example runners |
| `ultidock example run NAME [arguments]` | Invoke an example with its own options |
| `ultidock benchmark --help` | Discover benchmark subcommands |
| `ultidock report DIRECTORY` | Build report artifacts from recognized result files |
| `ultidock clean -- --help` | Inspect cleanup scope and opt-in artifact deletion |

## Pipeline options

`ultidock run -- --help` and `python docking/run.py --help` describe setup/docking
arguments. The separator forwards help to the underlying script; without it,
`--help` displays the CLI wrapper's own options. Wrapper commands such as `cavity`
accept pipeline arguments after
their own options. Common options include `--mode`, `--skip-wget`,
`--macro-mol-dir`, `--ligands-dir`, `--results-dir`, `--vina-cpu`,
`--vina-exhaustiveness`, `--vina-num-modes` and `--vina-seed`.
See [Configuration](configuration.md) for defaults and scope.

## Teacher and example options

```bash
ultidock example run quickstart -- --help
ultidock example run quickstart --dry-run
ultidock example run d2-antipsychotics -- --help
ultidock example run d2-antipsychotics p2rank --mode cpu
```

The teacher uses `--site-method`; the D2 runner uses an optional positional method.
`--yes` explicitly starts the teacher's real docking workflow without prompts.
The teacher and D2 `--dry-run` modes do not stage inputs or dock. Pipeline wrapper
dry runs can write configuration, and pocket-method dry runs can invoke predictors;
read each command's help rather than assuming all dry runs have identical effects.

## MolGuard

```bash
molguard doctor
molguard pdbqt check ligand.pdbqt
molguard pdbqt normalize ligand.pdbqt -o normalized.pdbqt
molguard receptor prepare source.pdb -o receptor.pdbqt
molguard grids --help
```

Receptor preparation and ligand normalization have different chemistry and record
handling. Choose the command for the file's actual role.
