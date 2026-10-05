# Changelog

## 1.1.2 - 2026-10-05

- Added the complete Read the Docs site using the official tutorial-template layout, with getting-started guides, scientific background, tutorials, reference and development documentation.
- Replaced the fixture-based quickstart with an interactive Next/Back/Exit teacher and a preview mode that does not run docking.
- Added the D2 receptor example with haloperidol, escitalopram and morphine, explicit input staging and CPU/GPU engine selection.
- Linked exported docking scores to the saved best-pose PDBQT, retaining input ligand identity and all parsed engine runs in SQLite.
- Improved receptor preparation and missing donor-hydrogen recovery, including heavy-atom preservation checks and preparation provenance.
- Renamed the Python distribution to `ultidock`, keeping the `molguard` command and imports.
- Included workflow modules, examples, guides and native build sources in wheels and source archives.
- Added `python -m ultidock` and a managed user workspace for regular installs, configurable with `ULTIDOCK_HOME`.
- Added installation checks outside the source checkout and preserved binding-site identifiers in analysis exports.
- Fixed manual known-site boxes being rejected when the receptor filename differs from the box label.
- Made pose-parsing tests independent of a locally generated docking configuration.

## 1.1.1 - 2026-10-01

- Added MolGuard checks and deterministic receptor preparation for molecular input files.
- Added known-site, CaV-EMPS cavity, blind, fpocket, and P2Rank docking workflows to the Ultidock CLI.
- Added a lightweight quickstart, run reports, and PyMOL/ChimeraX visualization helpers with command and version provenance.
- Added normalized COACH420/HOLO4K site-prediction benchmarking with DCA Top-n and Top-(n+2) metrics; targets with empty method output remain in the denominator.
- Improved CaV-EMPS site portfolios, adaptive cavity thresholds, grid box sizing, and receptor preparation diagnostics.
- Added parallel DUD-E cavity-recovery benchmarking and local pocket-tool installation support.
- Kept generated benchmark maps in temporary work directories by default to reduce disk use.
- Added the preregistered GABA-A 8DD2 case-study example and paper-facing COACH420 analysis documentation.
