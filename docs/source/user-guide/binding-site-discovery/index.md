---
orphan: true
---

# Binding-site discovery

A site finder proposes where to search. A docking engine then samples ligand
poses within each box. These are separate inference steps and require separate
validation. Choose the workflow using the information available before docking.

| Situation | Starting command |
| --- | --- |
| Trusted experimental/expert center | `ultidock known-site --center x,y,z` |
| Unknown site, receptor-based proposal | `ultidock cavity` |
| Alternative pocket methods | `ultidock fpocket` or `ultidock p2rank` |
| Broad whole-receptor baseline | `ultidock blind` |

The default input folders and automatic backend selection work for these
commands too. Pass directory overrides only when inputs are elsewhere. With
local ligands, add `--skip-wget` to avoid the bundled example download.
The number of sites is a search budget, not the number of true binding sites.
Retain all site identifiers and box dimensions with docking results. A site's
numeric label (S1, S2, …) is not a confidence statement.

```{toctree}
:maxdepth: 1

cav-emps
fpocket
p2rank
```
