# Binding-site prediction

A site predictor proposes regions to search. A docking engine then samples ligand
poses within those regions. These are distinct tasks: a well-localized site can
still yield an incorrect pose, and a favorable docking score cannot rescue a box
that excludes the experimental binding site.

Ultidock supports receptor-derived CaV-EMPS proposals and integrations with local
fpocket and P2Rank installations. A known-site run instead uses supplied coordinates.
A blind run uses a broad box; it is not equivalent to evaluating several local pockets.
See [Binding-site discovery](../user-guide/binding-site-discovery/index.md) for commands.

## Evaluate site recovery separately

For a reference complex, define the reference ligand atoms and site before running
the experiment. Measure center localization, ligand coverage or another explicit
criterion, and report the number of proposals allowed. Center-to-center distance
and atom coverage answer different questions. Report success at a stated cutoff
and top-k budget, rather than an unspecified “site accuracy.”

Withhold reference ligand coordinates from receptor-only prediction. Avoid tuning
thresholds on the same targets used for final reporting. Compare methods with the
same receptor preparation, box budget and evaluation definitions. A retained
cofactor, bound ligand or different chain selection can change the effective task.

The repository's site-prediction benchmarks keep prediction and evaluation stages
separate. Their manifests and evaluation settings are part of the result, not
optional bookkeeping.
