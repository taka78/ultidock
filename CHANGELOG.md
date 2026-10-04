# Changelog

## Unreleased

- Bound new AutoDock-GPU provenance to its coordinate DLG as well as its score XML, and recorded final MD preflight failures in the handoff status file.
- Added `--md-config`/`--md` continuation to all public docking modes, with MD preflight, selection restricted to the current docking run, automatic system preparation/equilibration and persistent handoff status. `ultidock md run --config ... --docking-dir ...` now also handles existing docking results in one command.
- Preserved high-throughput screening after individual ligand/receptor/site failures, with structured JSON/CSV failure reports and successful-output-only MD selection. Runs with no successes report failure without selecting stale poses. MD preflight failures disable MD while docking continues.
- Made SERT and 4COF examples default to automatic CPU fallback, accept MD continuation options, support isolated output directories and dry runs, and preserve caller-relative paths.
- Accelerated CaV-EMPS exterior flood filling while preserving its boundary seeds, voxel connectivity, cavity distances and component labels; verified exact equality on the bundled SERT grid.

- Unified docking, receptor preparation and MD dependency diagnostics in `ultidock doctor`; removed the separate `ultidock md doctor` command.

- Added CLI-only docking-to-GROMACS workflows under `md-simulation/`, including top-five distinct ligand selection, exact scored pose extraction, automatic GAFF2/AM1-BCC ligand parameters, soluble complexes and reviewed bilayer seed insertion.
- Added configurable water, box, salt, temperature and duration settings with force-field compatibility checks, restrained equilibration, checkpointed stages and explicit production review.
- Recorded docking output/input provenance and corrected new database `docking_file` values to reference the output pose container instead of the undocked input ligand.

- Renamed the Python distribution to `ultidock`, keeping the `molguard` command and imports.
- Included workflow modules, examples, guides and native build sources in wheels and source archives.
- Added `python -m ultidock` and a managed user workspace for regular installs, configurable with `ULTIDOCK_HOME`.
- Added installation checks outside the source checkout and preserved binding-site identifiers in analysis exports.
- Fixed manual known-site boxes being rejected when the receptor filename differs from the box label.

## 1.1.1 - 2026-10-01

- Added MolGuard checks and deterministic receptor preparation for molecular input files.
- Added known-site, CaV-EMPS cavity, blind, fpocket, and P2Rank docking workflows to the Ultidock CLI.
- Added a lightweight quickstart, run reports, and PyMOL/ChimeraX visualization helpers with command and version provenance.
- Added normalized COACH420/HOLO4K site-prediction benchmarking with DCA Top-n and Top-(n+2) metrics; targets with empty method output remain in the denominator.
- Improved CaV-EMPS site portfolios, adaptive cavity thresholds, grid box sizing, and receptor preparation diagnostics.
- Added parallel DUD-E cavity-recovery benchmarking and local pocket-tool installation support.
- Kept generated benchmark maps in temporary work directories by default to reduce disk use.
- Added the preregistered GABA-A 8DD2 case-study example and paper-facing COACH420 analysis documentation.
