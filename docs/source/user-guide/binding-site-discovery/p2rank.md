# P2Rank

P2Rank is an external machine-learning pocket predictor. Ultidock runs a local
installation and converts predictions into docking boxes. The installer targets
P2Rank 2.5; use a supported Java runtime (Java 21 is the documented setup choice).

```bash
bash scripts/install_pocket_tools.sh p2rank
ultidock doctor
ultidock example run sert-escitalopram p2rank --mode auto
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

P2Rank 2.5 needs Java 17–23; Java 25 can fail with
`Unsupported class file major version 69`. The bundled
`external/bin/prank` launcher selects Ubuntu Java 21 when `JAVA_HOME`
is unset. If `JAVA_HOME` points at an incompatible runtime, use:

```bash
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
"$JAVA_HOME/bin/java" -version
```

The default 35 Å box side and 0.375 Å grid spacing use the predictor's
reported center with no extra pocket padding. `sites.tsv` and raw
predictions stay in the run directory:

```bash
ultidock p2rank --receptor receptor.pdbqt --autosites 6 \
  --box-size 35
```

Missing P2Rank does not block setup for other modes; selecting it can
install the local pinned version. Pass `--receptor` explicitly when the
default receptor directory has multiple PDBQT files.
