# SERT–Escitalopram Docking Example

This example provides one prepared SERT receptor (`5i6x_edited.pdbqt`) and one
escitalopram ligand (`escitalopram-e.pdbqt`). Start with the installation and
GPU/CPU instructions in the [setup guide](../../SETUP.md). A full run needs
AutoGrid and either a working AutoDock-GPU backend or CPU-mode Vina.

From the repository root, run one site method at a time:

```bash
ultidock example run sert-escitalopram            # CaV-EMPS, the default
ultidock example run sert-escitalopram p2rank     # P2Rank 2.5
ultidock example run sert-escitalopram fpocket    # fpocket 4.2.3
```

With a source-checkout installation, replace `ultidock` with
`/usr/bin/python3 -m cli.ultidock`. If you have no GPU runtime, append
`--mode cpu` to any of the commands. The optional pocket programs are installed
locally when first selected; P2Rank needs Java 17–23, and fpocket needs the C
build dependencies in [requirements-ubuntu.txt](../../requirements-ubuntu.txt).

The runner prints a new directory under `workspace/<timestamp>/` for each
invocation. It stages only this example's receptor and ligand, skips the
general ligand download manifest, and keeps generated sites, grids, and docking
results in that directory. AutoGrid builds one map set for each proposed site,
so a multi-site run takes longer than pocket prediction alone. The local script
can also be run as `python3 example-run.py [cav-emps|p2rank|fpocket]` from this
folder.

Reference: [PubChem Compound Summary for Escitalopram](https://pubchem.ncbi.nlm.nih.gov/compound/Escitalopram).
