# P2Rank

P2Rank is an external machine-learning pocket predictor. Ultidock runs a local
installation and converts predictions into docking boxes. The installer targets
P2Rank 2.5; use a supported Java runtime (Java 21 is the documented setup choice).

```bash
bash scripts/install_pocket_tools.sh p2rank
ultidock doctor
ultidock example run d2-antipsychotics p2rank --mode cpu
```

To generate site boxes only:

```bash
ultidock pocket-box --method p2rank --receptor /absolute/receptor.pdbqt \
  --output /absolute/sites.tsv --autosites 6 --box-size 35
```

`--tool` overrides the executable. When `JAVA_HOME` selects an incompatible
runtime, fix the environment before interpreting a tool crash as a receptor
problem. The bundled launcher can choose an available supported Java when
`JAVA_HOME` is unset; check the actual runtime printed by the tool.

Retain the P2Rank predictions, version, box construction parameters and receptor
preparation provenance. Its confidence/ranking values are not interchangeable
with CaV-EMPS scores or Vina affinities. Combined `ultidock p2rank` mode accepts
pipeline arguments after the pocket options; `--dry-run` may run the predictor
before printing the proposed docking command.
