# Scoring / ranking

Docking scores are model-dependent estimates used to prioritize poses. They are
not measured affinities, and a difference between two values is not automatically
a meaningful biological effect. Compare like-for-like preparation, engines,
parameters and search budgets.

## Which pose is exported?

For Vina, Ultidock extracts the model with the lowest reported affinity into a
best-pose PDBQT. For AutoDock-GPU, the engine's saved `-best.pdbqt` can correspond
to a different run than the lowest binding-energy value in the XML, because the
engine's best-output selection uses its own energy criterion. Ultidock matches
saved coordinates to the DLG run before associating the pose with its score.

The analysis CSV's `docking_file` refers to the selected output PDBQT;
`ligand_file` refers to the prepared input. Check the associated run/model and site
when joining results. Other parsed poses remain in the database.

## Compare with care

- More sampled sites or poses give a ligand more opportunities to produce an extreme
  score. Keep the sampling budget consistent when comparing screens.
- Protonation, tautomer choice and receptor preparation can change scores. Retain
  these identities rather than collapsing distinct inputs under one compound name.
- Vina RMSD bounds relate its output poses to its reference/top pose; they do not
  establish agreement with an experimental crystal structure.
- Missing RMSD values in GPU results mean unavailable information, not a perfect match.
- Do not combine Vina scores, AutoDock energies and pocket-prediction scores into a
  single numerical ranking without a justified calibration procedure.

Inspect poses and use appropriate held-out labels or structural references to
assess enrichment or pose recovery. The three-ligand D2 example teaches the workflow;
it is not an enrichment benchmark. See [Results & reports](../user-guide/results-reports.md).
