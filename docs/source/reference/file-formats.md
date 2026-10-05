# File formats

| Artifact | Role |
| --- | --- |
| Receptor PDB / MOL2 | Structural source accepted by MolGuard conversion |
| Receptor PDBQT | Prepared coordinates, charges and AutoDock atom types |
| Ligand PDBQT | Prepared molecule with charges, types and ROOT/BRANCH/TORSDOF torsion records |
| `.pdbqt.gz` | Compressed molecular input; extraction behavior depends on workflow stage |
| `.prep.json` | Receptor preparation provenance and diagnostics |
| `.config.toml` | Optional per-receptor geometry profile |
| `centers.tsv` | Receptor/site identifiers, coordinates, dimensions and spacing |
| `.gpf`, `.glg`, `.fld`, `.map` | AutoGrid parameters, log, field descriptor and map data |
| `.xml`, `.dlg` | AutoDock-GPU run scores and docking log/coordinates |
| Output `.pdbqt`, `-best.pdbqt` | Engine poses and selected best pose |
| `ultidock_results.db` | SQLite docking records |
| `*-docking-results.csv` | Filtered analysis export |
| `run_config.yaml` | High-level command/run description |

## Site table

The current parser accepts whitespace-separated rows in this order:

```text
# receptor site_id cx cy cz nx ny nz spacing
receptor S1 12.3 4.5 -6.7 81 81 81 0.375
```

Coordinates and spacing are in Å; dimensions are integer grid counts. Receptor
keys normally match filename stems. Comment metadata controls policy/reuse and
manual-site handling, so preserve generated comments. Prefer `known-site` or
`pocket-box` to create site definitions rather than hand-editing dimensions.

## Input identity versus pose identity

An input `ligand.pdbqt` records preparation coordinates. An output PDBQT records
coordinates sampled by the engine, potentially in multiple MODEL blocks. CSV
`docking_file` should identify the selected output; `ligand_file` records the input.
Do not substitute one for the other because their filenames contain the same name.

PDBQT is not a lossless chemical exchange format. Retain the source SDF/MOL2 and
preparation metadata separately. See [Results & reports](../user-guide/results-reports.md)
for CSV and database interpretation.
