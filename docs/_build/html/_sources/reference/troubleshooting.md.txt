# Troubleshooting

Start with `ultidock doctor`, the printed workspace path and the first failing
subprocess log. Preserve its inputs and command before retrying.

| Symptom | Check / next action |
| --- | --- |
| `ultidock: command not found` | Activate the installation venv; from a checkout use `python -m cli.ultidock` |
| Teacher reports missing dependencies | Install listed Python/native tools using the installation guide, then restart; preview with `--dry-run` |
| Teacher requires a terminal | Run interactively, use `--dry-run` to read lessons, or explicitly use `--yes` to start docking in a script |
| No HD/HS donor hydrogens | Keep automatic receptor preparation enabled and install Open Babel; examine recovery warnings and `.prep.json` |
| Recovery fails heavy-atom preservation | Return to the source structure and prepare it with a suitable tool; do not bypass the check to force changed geometry through |
| Unknown atom type or non-finite coordinates/charges | Inspect the source molecule and preparation settings; normalization cannot reconstruct missing chemistry |
| GPU not detected | Check the allocated device, driver and CUDA/OpenCL runtime; choose `--mode cpu` if CPU execution is intended |
| AutoGrid fails | Read `.glg`; inspect receptor types, map parameters, executable and output permissions |
| P2Rank Java/class-file error | Select the supported Java runtime for the installed P2Rank version; inspect `JAVA_HOME` |
| Too many ligands discovered | Remove unrelated inputs from this run's ligand directory or use an isolated example runner |
| CSV is empty | Check docking failures and analysis filters; inspect unfiltered best-pose exports and the database |
| CSV points at an input ligand | Check whether rows are from an older database; new best-pose metadata is not retroactively reconstructed for all legacy records |
| HTML report has empty tables | Generic reporting requires recognized artifacts; inspect engine CSV/SQLite directly |
| Another run changes settings | Use separate application homes and output/grid directories for concurrent jobs |

## Receptor failures

MolGuard's automated recovery handles specified cases; it cannot make arbitrary
receptors chemically correct. Do not delete a preparation backup until you have
compared structures. A partially hydrogenated receptor may need full preparation
from its source even when the “no donor hydrogens” condition is absent.
See [Preparing a receptor](../user-guide/preparing-receptor.md).

## Useful bug reports

Include Ultidock and engine versions, command, active workspace, operating system,
backend, first relevant error and a minimal reproducible input when shareable.
Mention whether setup was skipped, whether a database was reused, and whether the
problem reproduces in a fresh example workspace. Avoid posting credentials or
restricted molecular datasets in public issues.
