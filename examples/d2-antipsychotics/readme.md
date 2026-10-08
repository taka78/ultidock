# D2 receptor: haloperidol and comparison ligands

Run the completed three-ligand example against `6CM4-edited.pdbqt`:

```bash
ultidock example run d2-antipsychotics --dry-run
ultidock example run d2-antipsychotics --mode auto
ultidock example run d2-antipsychotics --mode gpu
ultidock example run d2-antipsychotics fpocket --mode auto
ultidock example run d2-antipsychotics p2rank --mode auto
```

The D2 runner defaults to CPU Vina. `--mode auto` uses a detected GPU for
AutoDock-GPU, including distributing jobs across multiple NVIDIA GPUs, and
falls back to CPU Vina otherwise. CaV-EMPS is the default site finder.
`--mode gpu` requires a working GPU runtime; the other
site methods require their corresponding tools.
In a source checkout, replace `ultidock` with `/usr/bin/python3 -m cli.ultidock`.

| Input | Role |
| --- | --- |
| `6CM4-edited.pdbqt` | Edited D2 receptor input; the deposited structure contains risperidone. |
| `haloperidol.pdbqt` | Haloperidol, supplied from PubChem CID 3559. |
| `escitalopram-e.pdbqt` | Comparison ligand. |
| `morphine-e.pdbqt` | Comparison ligand. |

The runner stages exactly these files in a new `workspace/<timestamp>/`, prepares
the receptor, predicts sites, and docks all three ligands. Downloads are skipped.
Preparation may recover missing donor hydrogens in the receptor; inspect the
`.prep.json` report and any warnings before using the results scientifically.
Inputs in this folder are not overwritten by the run.

Inspect `RESULTS_DIR/` for SQLite and CSV results and `DOCKING_DIR/` for pose files.
The exported `docking_file` points to the saved best PDBQT pose with its matching
score. All runs remain in SQLite. These three molecules are a teaching/comparison
set, not a validated active/decoy benchmark; scores do not establish drug selectivity.

## Source notes supplied with the dataset

Structure of the D2 dopamine receptor bound to the atypical antipsychotic drug risperidone.
Wang, S., Che, T., Levit, A., Shoichet, B.K., Wacker, D., Roth, B.L.
(2018) Nature 555: 269-273

PubMed: 29466326 Search on PubMedSearch on PubMed Central
DOI: https://doi.org/10.1038/nature25758
Primary Citation of Related Structures:  
6CM4

PubMed Abstract: 
Dopamine is a neurotransmitter that has been implicated in processes as diverse as reward, addiction, control of coordinated movement, metabolism and hormonal secretion. Correspondingly, dysregulation of the dopaminergic system has been implicated in diseases such as schizophrenia, Parkinson's disease, depression, attention deficit hyperactivity disorder, and nausea and vomiting. The actions of dopamine are mediated by a family of five G-protein-coupled receptors. The D2 dopamine receptor (DRD2) is the primary target for both typical and atypical antipsychotic drugs, and for drugs used to treat Parkinson's disease. Unfortunately, many drugs that target DRD2 cause serious and potentially life-threatening side effects due to promiscuous activities against related receptors. Accordingly, a molecular understanding of the structure and function of DRD2 could provide a template for the design of safer and more effective medications. Here we report the crystal structure of DRD2 in complex with the widely prescribed atypical antipsychotic drug risperidone. The DRD2-risperidone structure reveals an unexpected mode of antipsychotic drug binding to dopamine receptors, and highlights structural determinants that are essential for the actions of risperidone and related drugs at DRD2.


National Center for Biotechnology Information. PubChem Compound Summary for CID 3559, Haloperidol. https://pubchem.ncbi.nlm.nih.gov/compound/Haloperidol. Accessed Oct. 5, 2026.
