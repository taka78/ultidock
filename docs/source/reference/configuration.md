# Configuration

The pipeline's setup stage writes `docking/config.py` in the active application
workspace. It is executable Python, not a user YAML schema. Wrapper commands also
write `run_config.yaml` as a run description; that file does not replace the
legacy pipeline configuration. Use `ultidock doctor` to find the active workspace.

## Public setup options

| Option | Default / meaning |
| --- | --- |
| `--mode` | `auto`; allows CPU fallback. `gpu` requires a detected GPU; `cpu`, `cuda`, `opencl` select explicitly |
| `--skip-wget` | Use local ligand inputs without executing the download list |
| `--macro-mol-dir`, `--ligands-dir` | Receptor and ligand directories under the application docking directory |
| `--docking-dir`, `--analysis-dir`, `--results-dir` | Output locations; use absolute paths for independent runs |
| `--vina-dir`, `--autodock-gpu-dir` | Native tool directories |
| `--grid-mode` | `centers`; also `ligand`, `residues`, `blind` |
| `--grid-spacing` | 0.375 Å |
| `--grid-margin` | 5.0 Å |
| `--grid-cap` | 150.0 Å |
| `--autosites` | 6 requested automatic sites |
| `--centers-tsv` | Site file path |
| `--ref-ligand-pdb` | Reference ligand for ligand-centered grids |
| `--vina-cpu` | 2 threads per Vina job |
| `--vina-exhaustiveness` | 8 |
| `--vina-num-modes` | 9 |
| `--vina-seed` | Unset unless supplied |
| `--receptor-prep-mode` | `auto`; `off` disables preparation |
| `--receptor-prepare-command` | Optional converter template with `{input}`, `{output}`, `{seed}` |
| `--receptor-prep-seed` | 42 |
| `--force-receptor-prep` | Regenerate same-stem receptor outputs |
| `--skip-profile` | Preserve existing receptor geometry sidecars |

High-level `known-site` accepts `--center` and `--box-size`; `cavity` accepts
`--autosites`. Their `--output-dir` groups output artifacts but does not move your
input receptor and ligand directories. Explicit input paths remain important.

## Advanced settings and concurrency

For a controlled setting without a CLI flag, first generate configuration using
`ultidock setup`, edit the active `docking/config.py`, then run
`ultidock run --skip-setup`. This preserves settings such as `SITE_POLICY`; setup
would otherwise overwrite them. With `--skip-setup`, do not expect new setup flags
to modify the existing configuration. Archive the actual configuration used.

`ULTIDOCK_HOME` selects a managed workspace. `ULTIDOCK_WORKERS` overrides worker
count; `GPU_SLOTS_PER_DEV` defaults to 2. Avoid sharing generated configuration or
writable receptor-grid directories between concurrent jobs. The engine subprocess
thread environment is set by the runner; see [HPC / batch screening](../user-guide/hpc-batch-screening.md).
