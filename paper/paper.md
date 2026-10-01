# Ultidock: a cavity-aware virtual screening workflow

## Statement of Need

Molecular docking workflows often require manual receptor preparation, manual
binding-site definition, fragile grid setup, and ad hoc post-processing. These
steps make large screening runs difficult to reproduce and difficult to compare
across automatic and known-site protocols.

Ultidock is a reproducible, cavity-aware virtual screening workflow that
automates receptor preparation, binding-site proposal, docking execution,
result storage, reporting, and visualization.

## Method Summary

Ultidock includes CaV-EMPS, Cavity detection via Electrostatic Map Pocket
Scoring, for receptor-only binding-site proposal. CaV-EMPS combines receptor
geometry, AutoGrid-derived map evidence, and calibrated site-priority ranking to
produce a compact set of docking boxes without using a co-crystallized ligand as
input.

CaV-EMPS ranking scores are used only to rank site hypotheses. They are not
binding affinities or docking scores.

## Benchmark Plan

Site prediction is evaluated with DCA-based Top-n and Top-(n+2) metrics on
COACH420. HOLO4K confirmation and DUD-E docking-enrichment evaluation remain
separate planned analyses.
using known-site, blind, external-pocket, and CaV-EMPS modes.

## Reproducibility

Ultidock run directories record configuration, command provenance, software
versions, generated site tables, reports, and visualization helper files.
