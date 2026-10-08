# GABAA Benzodiazepine Docking Example

This example stages the prepared 4COF receptor and three bundled ligands. Set
up Python, AutoGrid, and a docking backend using the [setup guide](../../SETUP.md)
before starting a full run.

From the repository root:

```bash
ultidock example run gabaa-benzos --mode auto
```

The runner defaults to GPU mode, which requires a working runtime. With
`--mode auto`, it uses AutoDock-GPU on a detected GPU and CPU Vina otherwise.
With multiple NVIDIA GPUs, jobs are distributed across the detected devices.
For a source-checkout installation, replace `ultidock` with
`/usr/bin/python3 -m cli.ultidock`. The runner creates a new directory under
`workspace/<timestamp>/`, copies only the bundled receptor and ligands, skips
the general ligand download manifest, and prints the directory path so you can
inspect the results. You can also run `python3 example-run.py` from this folder.

Structure source: [PDB 4COF](https://www.rcsb.org/structure/4COF).
