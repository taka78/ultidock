# What is Ultidock?

Ultidock is a molecular docking workflow for moving from prepared receptors and
ligands to binding-site hypotheses, docking poses, and result tables. It combines
Python orchestration with native scientific programs; installing its Python
package alone does not install a GPU driver or every native tool.

The main components are:

| Component | Responsibility |
| --- | --- |
| MolGuard | Check molecular file formatting, prepare receptors, and validate grids. |
| CaV-EMPS | Propose receptor-derived docking sites using geometry and AutoGrid maps. |
| fpocket / P2Rank adapters | Turn alternative pocket predictions into docking boxes. |
| AutoDock Vina / AutoDock-GPU | Search and score ligand poses. |
| SQLite and analysis tools | Preserve runs and export scores linked to pose files. |

Use a known-site box when you have an experimental or justified structural
reference. Use site discovery when the binding location is unknown. Whole-receptor
blind docking is another baseline, with a larger search space and different cost.

Ultidock is not an experimental affinity assay, a substitute for chemical model
curation, or proof of a drug's selectivity. Its output is a set of computational
hypotheses to inspect and validate. See [Limitations](../scientific-background/limitations.md).

Begin with the [interactive teacher](quick-start.md), or jump to the
[D2 example tutorial](../tutorials/small-molecule-docking.md).
