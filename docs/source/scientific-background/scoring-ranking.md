# Scoring and ranking

Docking scores are model-dependent estimates for prioritizing poses, rather
than measured affinities. Compare the same preparation, engine, parameters
and sampling budget before interpreting score differences.

## Which model does a row describe?

Ultidock records each parsed Vina MODEL or AutoDock-GPU XML run in SQLite.
Analysis exports models passing its filters, without automatic per-case
best-pose reduction. `docking_file` points to the coordinate container and
`model` identifies the scored pose within it. For AutoDock-GPU, retain both
the XML score file and DLG coordinates.

AutoDock-GPU's engine-selected `-best.pdbqt` can represent a different run from
the lowest XML binding-energy row. Use the row's run identity when inspecting
or comparing coordinates. Vina output retains scored MODEL blocks.

The optional MD handoff makes a separate selection: the lowest binding-energy
pose per distinct ligand across matching sites/models, restricted to the
protocol's receptor and engine. It ranks raw outputs independently of analysis
CSV filters. XML and DLG energy/run agreement is checked before GPU pose
preparation. See [Docking to MD](../user-guide/molecular-dynamics/docking-to-md.md).

## Compare with care

- More sites or poses offer more opportunities to obtain an extreme score.
  Keep sampling budgets consistent.
- Protonation, tautomers and receptor preparation can change scores; retain
  distinct input identities and preparation settings.
- Vina RMSD bounds compare its output poses with its reference/top pose,
  rather than an experimental crystal structure.
- Missing GPU RMSD values describe unavailable information.
- Vina scores, AutoDock energies and pocket-prediction scores have different
  meanings; combining them requires a justified calibration.

Inspect poses and use appropriate held-out labels or structural references
for enrichment or recovery. The three-ligand 4COF example demonstrates a
workflow, rather than an enrichment benchmark. MD trajectories also require
scientific interpretation; this branch does not automatically calculate binding
free energies. See [MD assumptions](molecular-dynamics.md) and
[Results and reports](../user-guide/results-reports.md).
