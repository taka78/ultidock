# Limitations

Ultidock automates preparation checks and engine orchestration. It cannot guarantee
that every receptor file represents a chemically valid or biologically relevant target.

## Structural and chemical assumptions

Missing residues or heavy atoms, ambiguous protonation, covalent ligands, unusual
cofactors and unsupported atom types require informed preparation. Existing
PDBQT canonicalization does not supply missing atoms or choose chemical states.
A file passing validation can still have incorrect chemistry. Docking uses the
supplied receptor structure and does not model arbitrary receptor flexibility,
membrane dynamics or solvent rearrangement.

## Molecular-dynamics scope

The optional MD workflow supports soluble systems and reviewed membrane seeds
with Amber-compatible protein force fields and GAFF2 ligands. It requires a
matched protein PDB, complete ligand SDFs and atom maps; docking inputs alone
cannot supply the system. Missing residues, unsupported chemistry and membrane
orientation/composition need review. Each protocol selects one receptor and one
scoring engine. MD complexes run sequentially; the docking multi-GPU worker
scheduler does not schedule GROMACS replicas.

New jobs stop at NPT by default, and production requires an equilibration review.
Completion of short test trajectories does not establish affinity, stability,
equilibration or a validated research protocol. See
[MD assumptions](molecular-dynamics.md) and [MD validation](../development/md-validation.md).

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
