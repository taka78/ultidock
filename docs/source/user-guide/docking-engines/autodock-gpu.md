# AutoDock-GPU

AutoDock-GPU uses receptor affinity maps and an installed CUDA or OpenCL runtime.
Build/runtime checks are part of setup; use explicit GPU mode when a CPU fallback
would invalidate your experiment:

```bash
ultidock example run d2-antipsychotics --mode gpu
```

The current runner requests 64 runs per ligand/site and enables XML and best-pose
output. These are runner settings, not a promise that `ultidock run` accepts every
native AutoDock-GPU flag. Check the command in logs when recording experiments.

| Output | Content |
| --- | --- |
| `.xml` | Run identifiers, energies and transformation data. |
| `.dlg` | Docked coordinates and detailed run/cluster information. |
| `-best.pdbqt` | The single best pose written by the engine. |

AutoDock-GPU selects its best file by total score; the first XML entry sorted by
binding energy can be a different run. Ultidock matches the best file's coordinates
to the DLG run before exporting its score. A missing or ambiguous match is reported
rather than linking the input ligand as if it were a docked pose.

GPU slots and worker concurrency should be chosen for available devices and
memory. See [HPC / batch screening](../hpc-batch-screening.md). Keep the native
engine version, parameter files and runtime information with any published result.
