# AutoDock Vina

Use CPU mode to run Vina:

```bash
ultidock example run sert-escitalopram --mode cpu
```

Vina receives the prepared receptor/ligand and a box derived from the selected
site. In the full Ultidock pipeline, AutoGrid remains part of receptor grid/site
preparation even when Vina performs the pose search.

The generated configuration includes `VINA_CPU`, `VINA_EXHAUSTIVENESS`,
`VINA_NUM_MODES`, and optional `VINA_SEED`. Set them on pipeline commands with
`--vina-cpu`, `--vina-exhaustiveness`, `--vina-num-modes`, and `--vina-seed`.
Defaults are 2 threads, exhaustiveness 8, up to 9 modes and no explicit seed.
These are setup options: with `--skip-setup`, the existing generated configuration
is reused. Keep the actual configuration with the run.

The output PDBQT contains MODEL blocks and `REMARK VINA RESULT` scores. Ultidock
stores parsed models in SQLite with the output container path and each model's
number and score. Analysis filters those rows without extracting a separate
best-pose PDBQT. Optional MD selection reads the scored MODEL coordinates.

Vina's RMSD lower/upper bounds describe differences from its reference/top pose
within the output, not an independent comparison to a crystal ligand. A negative
score is not proof of activity. See [Scoring / ranking](../../scientific-background/scoring-ranking.md)
and the [upstream Vina documentation](https://autodock-vina.readthedocs.io/).
