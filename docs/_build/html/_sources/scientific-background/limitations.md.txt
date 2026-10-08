# Limitations

Ultidock automates preparation checks and engine orchestration. It cannot guarantee
that every receptor file represents a chemically valid or biologically relevant target.

## Structural and chemical assumptions

Missing residues or heavy atoms, ambiguous protonation, covalent ligands, unusual
cofactors and unsupported atom types require informed preparation. Automated donor-
hydrogen recovery handles a specific input problem; a file passing validation can
still have incorrect chemistry. Standard runs use the supplied receptor structure
and do not model arbitrary receptor flexibility, membrane effects or solvent rearrangement.

## Sampling and scoring

Local boxes restrict the accessible search space. Broader boxes and additional
sites increase cost and may still miss relevant poses. Stochastic docking and
finite search budgets introduce variability. Scores should support prioritization,
not claims of efficacy, selectivity or experimental binding strength.

## Operational boundaries

The generated Python configuration is shared within an application workspace;
isolate concurrent jobs. GPU drivers and optional pocket predictors are external
dependencies. Large receptors can require substantial grid memory and disk space.
Generic HTML reports only populate tables from recognized input artifacts; consult
raw outputs and the database when a report is incomplete.

Record excluded targets, failed preparations, missing outputs and empty filtered
exports. Treat failures as part of a benchmark's outcome rather than silently
removing them. For actionable errors, consult [Troubleshooting](../reference/troubleshooting.md).
