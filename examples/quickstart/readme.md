# Ultidock Quickstart

This is the recommended first run:

Install the Python dependencies and CLI in [SETUP.md](../../SETUP.md#step-1--install-prerequisites)
first. From the repository root:

```bash
ultidock example run quickstart
```

For a source-checkout installation, use
`/usr/bin/python3 -m cli.ultidock example run quickstart`. This command does
not need Java, fpocket, P2Rank, AutoGrid, a GPU, or a vendor compute runtime.

The quickstart is intentionally lightweight. It creates the standard Ultidock
researcher-facing artifact set without requiring a full docking toolchain on the
first command:

- `input/receptor.pdb`
- `input/reference_ligand.mol2`
- `run_config.yaml`
- `sites.tsv`
- `predictions.tsv`
- `top_hits.csv`
- `results.sqlite`
- `report.md`
- `report.html`
- PyMOL and ChimeraX helper files

The site method shown in the report is `cav-emps`, the CaV-EMPS method:
Cavity detection via Electrostatic Map Pocket Scoring.

The bundled receptor and ligand are tiny synthetic public-domain fixtures used
for artifact and report smoke testing. The ligand is recorded as reference
metadata only; it is not used as an input to CaV-EMPS site generation.
