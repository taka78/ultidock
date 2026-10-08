---
orphan: true
---

# Docking engines

Ultidock's `--mode` selects a computational backend. CPU uses AutoDock Vina;
CUDA/OpenCL use AutoDock-GPU. They are distinct engines and scoring/search paths,
so do not pool their raw scores as though they were identical measurements.

| Mode | Behavior |
| --- | --- |
| `cpu` | Select Vina explicitly. |
| `gpu` | Require a detected supported GPU. |
| `cuda` | Request the CUDA backend. |
| `opencl` | Request the OpenCL backend. |
| `auto` | Permit detected GPU selection or CPU fallback. |

A successfully installed Python CLI does not establish that native programs or
device runtimes work. Run `ultidock doctor` and a small example before scaling.

```{toctree}
:maxdepth: 1

vina
autodock-gpu
```
