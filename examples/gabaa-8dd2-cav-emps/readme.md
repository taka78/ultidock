# GABAA 8DD2 Blind Site-Recovery Case Study

This example runs the preregistered case study **Blind recovery of
pharmacologically distinct GABAA-receptor binding regions** on PDB 8DD2. It is
separate from the older `gabaa-benzos` docking example, which uses 4COF and a
different scientific question.

## Scientific Question

Can a frozen six-box CaV-EMPS portfolio recover five experimentally occupied,
chemically distinct regions of the human alpha1-beta2-gamma2 GABAA receptor
without seeing any ligand coordinates?

The 2.90 Angstrom 8DD2 model contains two orthosteric GABA molecules, one
extracellular alpha1/gamma2 zolpidem molecule, and two beta2/alpha1
transmembrane zolpidem molecules. Prediction input retains only receptor chains
A-E. The script excludes all HETATM records, glycans, and Fab chains I-L before
any method runs. Ground-truth ligands are extracted into a separate directory
that is never passed to CaV-EMPS, fpocket, or P2Rank.

The checked-in [`preregistration.json`](preregistration.json) freezes the five
ligand instances, seven CaV-EMPS profiles, top-six portfolio, 35.625 Angstrom
standardized box, containment buffers, and physical hypotheses before results
are inspected. The seven profiles are the repository's controlled suite:
`combined-final` plus six one-component controls.

## Run

From the repository root:

```bash
ultidock example run gabaa-8dd2-cav-emps \
  --output-dir /path/to/gabaa_8dd2_case \
  --jobs 1 \
  --force
```

The script downloads the exact RCSB PDB revision recorded in the
preregistration file. For an offline or archived input:

```bash
ultidock example run gabaa-8dd2-cav-emps \
  --pdb /path/to/8DD2.pdb \
  --output-dir /path/to/gabaa_8dd2_case \
  --force
```

External methods can be located explicitly:

```bash
ultidock example run gabaa-8dd2-cav-emps \
  --output-dir /path/to/gabaa_8dd2_case \
  --fpocket-cmd /path/to/fpocket \
  --p2rank-cmd "/path/to/prank predict" \
  --force
```

Useful staged operations are `--prepare-only` and `--evaluate-only`. A normal
run executes preparation, all selected predictors, evaluation, and reporting.
Use `--methods cav-emps` for the CaV-EMPS profiles alone.

## Primary Endpoint

The primary endpoint is full-ligand containment by an axis-aligned docking box,
not center distance. A ligand is contained with buffer `b` only when every
ligand atom lies at least `b` Angstrom inside every box face. Results are
reported at 0, 2, and 4 Angstrom buffers for:

- native CaV-EMPS boxes;
- identical 35.625 Angstrom boxes around every method's top-six centers.

Distance to the closest ligand atom and distance to the ligand centroid are
secondary diagnostics. Raw C/e/d values sampled at a predicted center are
descriptive map values, not isolated energetic contributions. The controlled
ablations provide the stronger evidence about channel contribution.

## Outputs

```text
<output-dir>/
  source/8DD2.pdb
  normalized/gabaa_8dd2/8dd2/receptor_input.pdb
  ground_truth/labels.json
  ground_truth/ligands/*.pdb
  predictions/{cav-emps,fpocket,p2rank}/...
  evaluation/ground_truth_sites.csv
  evaluation/predictions_all.tsv
  evaluation/per_box_containment.csv
  evaluation/per_site_interpretation.csv
  evaluation/portfolio_summary.csv
  evaluation/ablation_summary.csv
  reports/case_study_report.md
  reports/case_study_report.html
  reports/containment_matrix.svg
  visualization/gabaa_8dd2_case_study.pml
  provenance.json
```

Missing tools or empty method output remain explicit in the summary tables;
they are never silently removed from the denominator.

## Interpretation Limits

The C/e/d expectations in the preregistration file are hypotheses made from
pocket chemistry. The cited structural literature establishes ligand locations
and chemistry, while the controlled ablations test what each CaV-EMPS channel
actually contributes.

The receptor input has no explicit lipid bilayer. Several GABAA transmembrane
pockets are partly lipid-facing, so this is both a difficult stress test and an
honest limitation. The reported additional gamma2/beta2 diazepam TMD site can
help interpret a sixth prediction, but it is not one of the five zolpidem/GABA
ground-truth labels. This case study does not replace the Holo4K containment
audit.

## Structural References

- [RCSB PDB 8DD2](https://www.rcsb.org/structure/8DD2), PDB DOI
  `10.2210/pdb8DD2/pdb`.
- Masiulis et al. (2019), *GABAA receptor signalling mechanisms revealed by
  structural pharmacology*, DOI `10.1038/s41586-018-0832-5`.
- Kim et al. (2020), *Shared structural mechanisms of general anaesthetics and
  benzodiazepines*, DOI `10.1038/s41586-020-2654-5`.
- Zhu et al. (2022), *Structural and dynamic mechanisms of GABAA receptor
  modulators with opposing activities*, DOI `10.1038/s41467-022-32212-4`.
- Legesse et al. (2023), *Structural insights into opposing actions of
  neurosteroids on GABAA receptors*, DOI `10.1038/s41467-023-40800-1`.
