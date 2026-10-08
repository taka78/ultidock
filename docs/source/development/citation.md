# Citation

Use the repository's [CITATION.cff](https://github.com/taka78/ultidock/blob/gmx-dev/CITATION.cff)
as the source for software citation metadata. Cite the release and commit actually
used, together with the workflow parameters needed to reproduce the experiment.

The current repository metadata identifies:

> Turgut, Taha. Ultidock, version 1.1.1 (2026-10-01). Software.
> https://github.com/taka78/ultidock

The citation file does not currently supply a DOI for an Ultidock methods paper;
do not invent one. A software citation does not replace citations for the docking
engine, site predictor, preparation software or structural dataset you used.
Consult those projects' own citation guidance for the relevant versions.

This branch adds MD while retaining the version recorded in `CITATION.cff`.
Record the `gmx-dev` commit used so the simulation code is identifiable.
For examples, retain the input provenance in each directory's README and
`dataset.json` when provided. For MD, cite the GROMACS version, AmberTools/ACPYPE,
protein and lipid force fields, water model and source structures used. Record
system preparation, chemical states, membrane assumptions and protocol settings
with the [MD artifacts](../user-guide/molecular-dynamics/running-and-results.md).

Ultidock is released under the
[MIT License](https://github.com/taka78/ultidock/blob/gmx-dev/LICENSE).
When applicable, also cite the engine methods used:

- Trott, O., and Olson, A. J. (2010). *AutoDock Vina: Improving the speed and
  accuracy of docking with a new scoring function, efficient optimization,
  and multithreading*. Journal of Computational Chemistry 31(2), 455–461.
- Santos-Martins, D., et al. (2021). *Accelerating AutoDock4 with GPUs and
  Gradient-Based Local Search*. Journal of Chemical Theory and Computation
  17(2), 1060–1073.
