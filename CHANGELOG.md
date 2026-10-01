# Changelog

## 1.1.1 - 2026-10-01

- Added MolGuard checks and deterministic receptor preparation for molecular input files.
- Added known-site, CaV-EMPS cavity, blind, fpocket, and P2Rank docking workflows to the Ultidock CLI.
- Added a lightweight quickstart, run reports, and PyMOL/ChimeraX visualization helpers with command and version provenance.
- Added normalized COACH420/HOLO4K site-prediction benchmarking with DCA Top-n and Top-(n+2) metrics; targets with empty method output remain in the denominator.
- Improved CaV-EMPS site portfolios, adaptive cavity thresholds, grid box sizing, and receptor preparation diagnostics.
- Added parallel DUD-E cavity-recovery benchmarking and local pocket-tool installation support.
- Kept generated benchmark maps in temporary work directories by default to reduce disk use.
- Added the preregistered GABA-A 8DD2 case-study example and paper-facing COACH420 analysis documentation.
