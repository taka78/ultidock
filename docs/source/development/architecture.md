# Architecture

Ultidock separates user commands, workspace management, molecular validation,
engine orchestration and result interpretation.

| Area | Responsibility |
| --- | --- |
| `cli/ultidock.py` | High-level modes, examples, benchmarks and reports |
| `cli/molguard.py` | Molecular validation/preparation commands |
| `ultidock/paths.py`, `ultidock/runtime_files.py` | Active workspace and packaged runtime resources |
| `docking/run.py`, `docking/setup.py` | Pipeline entry and generated configuration |
| `docking/dock_v02.py` | Receptor/site loading, worker dispatch and engine output handling |
| `docking/make_grids.py` | Site generation and grid preparation |
| `docking/pocket_boxes.py` | Local fpocket/P2Rank prediction adapters |
| `molguard/io/` | Fixed-format I/O, receptor preparation and validation |
| `docking/result_poses.py` | Output-pose extraction and score/coordinate association |
| `docking/db_manager.py` | SQLite persistence |
| `docking/analyse_docking_results.py` | Filtered analysis exports |
| `examples/common.py` | Isolated example staging and pipeline invocation |
| `benchmarks/` | Dataset, prediction and docking evaluation workflows |

## Execution and persistence

The full pipeline performs setup, extraction, docking and analysis. Setup writes a
Python configuration consumed by legacy pipeline modules. Receptor preparation
precedes grid construction; worker jobs then combine receptor/site/ligand inputs.
Keep configuration mutations outside parallel job execution.

An installed package materializes scripts, example data and native build sources
in a managed writable workspace, rather than writing outputs into site-packages.
Source-checkout behavior and `ULTIDOCK_HOME` selection are covered by packaging tests.

## Documentation build

The Read the Docs tutorial-template layout keeps Sphinx sources in `docs/source/`
and generated output in `docs/build/`. The root `index.rst` supplies navigation;
Sphinx reads the Markdown guides through MyST without importing the scientific
runtime. `.readthedocs.yaml` installs only `docs/requirements.txt` and selects
`docs/source/conf.py`. Add new pages to the appropriate toctree in
`docs/source/index.rst` or its nested guide indexes. `docs/Makefile` and
`docs/make.bat` provide the template's local build commands.
