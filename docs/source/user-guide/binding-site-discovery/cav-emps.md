# CaV-EMPS

CaV-EMPS (Cavity detection via Electrostatic Map Pocket Scoring) combines receptor
geometry with receptor-derived AutoGrid interaction maps. It is the default site
method in the D2 and SERT example runners.

```bash
ultidock example run d2-antipsychotics --mode auto
```

For your own prepared directories:

```bash
ultidock cavity --autosites 6 \
  --macro-mol-dir /absolute/inputs/receptors \
  --ligands-dir /absolute/inputs/ligands \
  --output-dir /absolute/runs/cav-emps --skip-wget
```

The pipeline builds receptor maps, constructs internal and surface candidates,
and selects separated boxes under the configured site policy. Keep `sites.tsv`,
AutoGrid `.gpf`/`.glg` files and per-receptor overrides with the result.

`R_MIN_CAVITY_A=None` enables adaptive cavity-radius selection. `AUTOSITES`,
`HOTSPOT_BOX_ANGLE`, separation bounds and surface-shell settings affect coverage
and cost. Record changes rather than tuning against an evaluation ligand without
reporting that information use.

CaV-EMPS site scores rank candidates; they are not docking affinities. Consult
[Methodology](../../scientific-background/cav-emps-methodology.md) for the signals
and [Limitations](../../scientific-background/limitations.md) before interpreting
site-recovery or enrichment measurements.
