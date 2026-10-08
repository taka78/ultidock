# Docking workflow

A run is a chain of decisions: molecular preparation → search-region selection →
pose sampling → scoring → interpretation. A failure early in that chain can make
plausible-looking scores misleading.

1. **Prepare molecular states.** Select the receptor assembly and chemical states
   of the ligands. Preserve source structures and preparation reports. MolGuard
   validates engine inputs; validation does not establish biological suitability.
2. **Define search regions.** Supply a known site, ask a pocket method for candidate
   sites, or use a whole-receptor box. The choice determines which poses are reachable.
3. **Prepare maps and sample poses.** Ultidock builds the required AutoGrid resources
   and dispatches receptor/site/ligand jobs to the chosen engine. Each engine has
   its own sampling and scoring model.
4. **Associate coordinates and scores.** Keep the actual output pose attached to its
   engine run/model. Ultidock retains raw outputs and stores parsed records in SQLite.
5. **Inspect and evaluate.** Review preparation warnings, placements and interactions.
   Use independent reference structures or experimental labels when evaluating accuracy.

For reproducibility, archive the source inputs, prepared files, site definitions,
configuration, engine versions, random seeds, logs and raw outputs. A fixed seed
helps repeat a specified setup; it does not establish cross-engine equivalence.
See [Scoring / ranking](scoring-ranking.md) before interpreting a hit list.
