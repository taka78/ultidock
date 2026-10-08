# File formats

| Artifact | Role |
| --- | --- |
| Receptor PDB / MOL2 | Structural source accepted by receptor conversion |
| Receptor PDBQT | Prepared coordinates, charges and AutoDock atom types |
| Ligand PDBQT | Prepared molecule with charges, types and torsion records |
| `.pdbqt.gz` | Compressed molecular input; extraction depends on the workflow stage |
| `.config.toml` | Optional per-receptor geometry profile |
| `centers.tsv` | Receptor/site IDs, coordinates, dimensions and spacing |
| `.gpf`, `.glg`, `.fld`, `.map` | AutoGrid parameters, log, field descriptor and maps |
| Vina output `.pdbqt` | Scored coordinates in MODEL blocks |
| AutoDock-GPU `.xml`, `.dlg` | Paired run scores and coordinates |
| AutoDock-GPU `-best.pdbqt` | Engine-selected single pose; does not identify every scored run |
| `ultidock_results.db` | SQLite docking records |
| `*-docking-results.csv` | Filtered model analysis |
| `run_config.yaml` | High-level command/run description |
| `docking-run-*.json`, `*.failures.csv` | Current successful-output manifest and case-failure report |
| Pose `.md.json` sidecar | Hashed receptor/ligand/output provenance for MD selection |
| `docking-run-*.md.json` | Integrated MD handoff status and job path |
| MD protocol `.json` | Reviewed chemistry, force field, system and execution settings |
| Ligand `.sdf` plus atom map | Chemical identity linked to selected docking coordinates |
| MD `.gro`, `.top`, `.itp` | Coordinates and system/molecule topology |
| MD `.mdp`, `.tpr` | Stage parameters and compiled simulation input |
| MD `.cpt`, `.xtc`, `.edr` | Checkpoint, compressed trajectory and energies |

## Site table

The parser accepts whitespace-separated rows in this order:

```text
# receptor site_id cx cy cz nx ny nz spacing
receptor S1 12.3 4.5 -6.7 81 81 81 0.375
```

Coordinates and spacing are in Å; dimensions are integer grid counts. Receptor
keys normally match filename stems. Generated comment metadata controls
policy/reuse and manual-site handling. Prefer `known-site` or `pocket-box` to
create definitions rather than hand-editing dimensions.

## Input, pose and MD identity

An input ligand PDBQT records preparation coordinates. A Vina output can contain
multiple models; an AutoDock-GPU DLG contains multiple runs. `docking_file` and
`model` together identify a scored coordinate set in this branch.

PDBQT is not a lossless chemical exchange format. Preserve source SDF/MOL2 and
chemical preparation separately. MD atom maps link PDBQT serial numbers to
zero-based SDF atom indices; membrane rigid transforms use Å while GRO
coordinates and membrane bounds use nm. See [MD protocol](md-protocol.md),
[Results and reports](../user-guide/results-reports.md) and
[MD artifacts](../user-guide/molecular-dynamics/running-and-results.md).
