# Why Ultidock?

Ultidock exists for researchers who need a reproducible path from receptor and
ligand inputs to docking results, without hand-editing grid boxes, receptor
formats, and post-run tables for every target.

The central automation problem is binding-site definition. Many screening tasks
do not start with a co-crystallized ligand or a trusted manual box. Ultidock's
CaV-EMPS method, Cavity detection via Electrostatic Map Pocket Scoring, proposes
receptor-derived docking sites from geometry and AutoGrid map evidence so that
automatic workflows can be compared against known-site and blind baselines.

Ultidock is designed around four claims:

- Deterministic receptor handling should be shared by normal runs and benchmarks.
- Automatic site proposal should be evaluated separately from docking enrichment.
- Reports and visual exports should be first-class outputs, not afterthoughts.
- Every benchmark should leave enough provenance to rerun or criticize it.

CaV-EMPS ranking scores are site-priority scores. They are not binding energies,
affinities, or pose scores.
